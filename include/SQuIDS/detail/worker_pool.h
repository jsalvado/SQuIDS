#ifndef SQUIDS_DETAIL_WORKER_POOL_H
#define SQUIDS_DETAIL_WORKER_POOL_H

// A small fork-join pool for splitting a loop over the nodes of the right-hand
// side across threads. The calling thread does one share of the work itself.
//
// Right-hand-side evaluations come microseconds apart during an integration,
// faster than a condition variable can wake a thread, so idle workers spin for a
// short while before going to sleep; once an integration ends they sleep and
// cost nothing.

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace squids{
namespace detail{

class WorkerPool {
public:
  typedef std::function<void(unsigned int, unsigned int, unsigned int)> Job; // (begin, end, thread)

  ///\param n total number of threads, including the caller
  explicit WorkerPool(unsigned int n): nthreads(n < 1 ? 1 : n){
    for(unsigned int id = 1; id < nthreads; id++)
      threads.emplace_back(&WorkerPool::Worker, this, id);
  }

  ~WorkerPool(){
    {
      std::lock_guard<std::mutex> lock(m);
      quit.store(true, std::memory_order_release);
      epoch.fetch_add(1, std::memory_order_release);
    }
    cv.notify_all();
    for(std::thread& t : threads)
      t.join();
  }

  WorkerPool(const WorkerPool&) = delete;
  WorkerPool& operator=(const WorkerPool&) = delete;

  unsigned int Size() const { return nthreads; }

  ///\brief Calls job(begin, end, thread) on contiguous pieces of [0, n) and
  /// returns when all are done. Threads claim small chunks dynamically, so fast
  /// and slow cores (e.g. performance and efficiency cores) share the work
  /// evenly. Rethrows the first exception.
  void Run(unsigned int n, const Job& job){
    if(nthreads == 1 || n < 2){
      job(0, n, 0);
      return;
    }
    std::lock_guard<std::mutex> one_job_at_a_time(run_m); // pools are shared
    current = &job;
    items = n;
    chunk = std::max(1u, n/(8*nthreads));
    next.store(0, std::memory_order_relaxed);
    error = nullptr;
    pending.store(nthreads-1, std::memory_order_relaxed);
    {
      std::lock_guard<std::mutex> lock(m);
      epoch.fetch_add(1, std::memory_order_release);
    }
    cv.notify_all();
    RunShare(0);
    while(pending.load(std::memory_order_acquire) != 0)
      std::this_thread::yield();
    if(error)
      std::rethrow_exception(error);
  }

private:
  unsigned int nthreads;
  std::vector<std::thread> threads;
  std::atomic<unsigned long> epoch{0};
  std::atomic<unsigned int> pending{0};
  std::atomic<bool> quit{false};
  const Job* current = nullptr;
  unsigned int items = 0;
  unsigned int chunk = 1;
  std::atomic<unsigned int> next{0};
  std::exception_ptr error;
  std::mutex m, error_m, run_m;
  std::condition_variable cv;

  void RunShare(unsigned int id){
    try{
      for(;;){
        const unsigned int begin = next.fetch_add(chunk, std::memory_order_relaxed);
        if(begin >= items)
          break;
        (*current)(begin, std::min(items, begin+chunk), id);
      }
    } catch(...){
      std::lock_guard<std::mutex> lock(error_m);
      if(!error)
        error = std::current_exception();
      next.store(items, std::memory_order_relaxed); // stop handing out work
    }
  }

  void Worker(unsigned int id){
    unsigned long seen = 0;
    for(;;){
      // spin briefly: the next evaluation usually follows within microseconds
      unsigned long e = epoch.load(std::memory_order_acquire);
      for(int spin = 0; e == seen && spin < 4000; spin++){
        std::this_thread::yield();
        e = epoch.load(std::memory_order_acquire);
      }
      if(e == seen){
        std::unique_lock<std::mutex> lock(m);
        cv.wait(lock, [&]{ return epoch.load(std::memory_order_acquire) != seen; });
        e = epoch.load(std::memory_order_acquire);
      }
      seen = e;
      if(quit.load(std::memory_order_acquire))
        return;
      RunShare(id);
      pending.fetch_sub(1, std::memory_order_acq_rel);
    }
  }
};

///\brief The process-wide pool with n threads, created on first use and torn
/// down when nothing holds it. Objects asking for the same number of threads
/// share it, so many SQuIDS objects (e.g. the zeniths of nuSQUIDSAtm) do not
/// each start their own threads.
inline std::shared_ptr<WorkerPool> SharedWorkerPool(unsigned int n){
  static std::mutex m;
  static std::map<unsigned int, std::weak_ptr<WorkerPool> > pools;
  std::lock_guard<std::mutex> lock(m);
  std::shared_ptr<WorkerPool> p = pools[n].lock();
  if(!p){
    p = std::make_shared<WorkerPool>(n);
    pools[n] = p;
  }
  return p;
}

} // namespace detail
} // namespace squids

#endif

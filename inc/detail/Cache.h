#ifndef SQUIDS_CACHE_H
#define SQUIDS_CACHE_H

#include <atomic>

namespace squids{
namespace detail{
 
///A simple, lock-free, thread safe object cache
///Heavily inspired by Chris Wellons' lock free stack:
/// http://nullprogram.com/blog/2014/09/02/
///\tparam T The type of object being stored. Must be default constructable, and
///          copyable. Default constructed instances will be used to communicate
///          lack of an object.
///\tparam N The maximum number of entries to store in the cache.
///\todo This class is now messy, since in our use case we don't actually need it
///      to be shared, so we can give each thread its own cache if the compiler
///      supports it. In that case, the cache can be much simpler (and somewhat
///      faster), since it knows it will only ever be accessed by one thread, so
///      no synchronization is required. This should probably be split into two
///      separate classes at some point, a simple cache and a shared cache. 
template<typename T, unsigned int N>
class cache{
private:
  ///One object record within a linked list
  struct record{
    record* next;
    T data;
    operator T(){
      return(data);
    }
  };
  
  static constexpr size_t max_buffer_size=N;
  record entries[max_buffer_size];
  
  ///The head of a linked list.
  ///Must be compact so that it can be updated atomically.
  struct list_head{
#ifndef SQUIDS_THREAD_LOCAL
    uint32_t counter; //data structure version number
#endif
    uint32_t index; //index of head node within entries, max_buffer_size indicates none
  };
#ifdef SQUIDS_THREAD_LOCAL
  using ListType=list_head;
#else
  using ListType=std::atomic<list_head>;
#endif
  ListType free_list, data_list;
  
  ///Remove the first entry from a linked list, if the list is not empty
  ///\param list The list to modify
  ///\return The address of the popped entry, or nullptr if the list was empty
  record* pop(ListType& list){
#ifdef SQUIDS_THREAD_LOCAL
    list_head orig=list;
    if(list.index==max_buffer_size) //empty stack
      return(nullptr);
    auto next_ptr=entries[list.index].next;
    if(next_ptr)
      list.index=next_ptr-entries;
    else
      list.index=max_buffer_size;
    return(entries+orig.index);
#else
    list_head orig=list.load(), next;
    do{
      if(orig.index==max_buffer_size) //empty stack
        return(nullptr);
      next.counter=orig.counter+1;
      auto next_ptr=entries[orig.index].next;
      if(next_ptr)
        next.index=next_ptr-entries;
      else
        next.index=max_buffer_size;
    }while(!std::atomic_compare_exchange_weak(&list, &orig, next));
    return(entries+orig.index);
#endif
  }
  
  ///Append an entry to a linked list
  ///\param list The list to modify
  ///\param The address of the entry to append
  void push(ListType& list, record* node){
#ifdef SQUIDS_THREAD_LOCAL
    if(list.index==max_buffer_size)
      node->next=nullptr;
    else
      node->next=entries+list.index;
    list.index=node-entries;
#else
    list_head orig=list.load(), next;
    uint32_t idx=node-entries;
    next.index=idx;
    do{
      next.counter=orig.counter+1;
      if(orig.index==max_buffer_size)
        node->next=nullptr;
      else
        node->next=entries+orig.index;
    }while(!std::atomic_compare_exchange_weak(&list, &orig, next));
#endif
  }
  
public:
  cache(){
    //initially, all entries are in the free list
    record* ptr=nullptr;
    for(unsigned int i=0; i<max_buffer_size; i++){
      entries[i].next=ptr;
      ptr=&entries[i];
    }
#ifdef SQUIDS_THREAD_LOCAL
    free_list=list_head{max_buffer_size-1};
    //no entries are in the data list
    data_list=list_head{max_buffer_size};
#else
    free_list.store({0,max_buffer_size-1});
    //no entries are in the data list
    data_list.store({0,max_buffer_size});
#endif
  }
  
  ///Insert an object into the cache if it is not full.
  ///\param value The object to store
  ///\return true if the entry was successfully inserted, or false if
  ///        insertion failed and the storage should be deallocated by
  ///        the caller
  bool insert(T value){
    record* entry=pop(free_list);
    if(!entry){ //no space left in cache
      return(false);
    }
    entry->data=value;
    push(data_list,entry);
    return(true);
  }
  
  ///Fetch an object from the cache if it is not empty.
  ///\return an available pointer/offset pair if possible, or a pair
  ///        {nullptr,0} if no storage is available in the cache
  T get(){
    record* entry=pop(data_list);
    if(!entry) //no cached memory available
      return(T());
    push(free_list,entry);
    return(*entry);
  }
};
  
}
} //namespace squids

#endif //SQUIDS_CACHE_H
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
    uint32_t counter; //data structure version number
    uint32_t index; //index of head node within entries, max_buffer_size indicates none
  };
  std::atomic<list_head> free_list;
  std::atomic<list_head> data_list;
  
  ///Remove the first entry from a linked list, if the list is not empty
  ///\param list The list to modify
  ///\return The address of the popped entry, or nullptr if the list was empty
  record* pop(std::atomic<list_head>& list){
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
  }
  
  ///Append an entry to a linked list
  ///\param list The list to modify
  ///\param The address of the entry to append
  void push(std::atomic<list_head>& list, record* node){
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
  }
  
public:
  cache(){
    //initially, all entries are in the free list
    record* ptr=nullptr;
    for(unsigned int i=0; i<max_buffer_size; i++){
      entries[i].next=ptr;
      ptr=&entries[i];
    }
    free_list.store({0,max_buffer_size-1});
    //no entries are in the data list
    data_list.store({0,max_buffer_size});
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
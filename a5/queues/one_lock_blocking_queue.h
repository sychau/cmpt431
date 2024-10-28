#include "../common/allocator.h"

template <class T>
class Node
{
};

extern std::atomic<bool> no_more_enqueues;

template <class T>
class OneLockBlockingQueue
{
    CustomAllocator my_allocator_;
public:
    OneLockBlockingQueue() : my_allocator_()
    {
        std::cout << "Using OneLockBlockingQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        //my_allocator_.freeNode(newNode);
    }

    void enqueue(T value)
    {
	    // add your enqueue code here
    }

    bool dequeue(T *value)
    {
        bool ret_value = false;
        return ret_value;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }
};

#include "../common/allocator.h"
#include <mutex>

template <class T>
class Node
{
public:
    T value;
    Node<T>* next;
};

template <class T>
class TwoLockQueue
{
    Node<T>* q_head_; // First node is sentinel
    Node<T>* q_tail_;
    CustomAllocator my_allocator_;
    std::mutex dequeue_m;
    std::mutex enqueue_m;

public:
    TwoLockQueue() : q_head_(nullptr), q_tail_(nullptr), my_allocator_() 
    {
        std::cout << "Using TwoLockQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));

        // Initialize the queue head or tail here
        Node<T>* sentinel = (Node<T>*)my_allocator_.newNode();
        sentinel->next = nullptr;
        q_head_ = sentinel;
        q_tail_ = sentinel;
    }

    void enqueue(T value)
    {
        enqueue_m.lock();
        Node<T>* new_node = (Node<T>*)my_allocator_.newNode();
        new_node->value = value;
        new_node->next = nullptr;

        q_tail_->next = new_node;
        q_tail_ = new_node;
        enqueue_m.unlock();
    }

    bool dequeue(T *value)
    {
        dequeue_m.lock();
        Node<T>* node = q_head_;
        Node<T>* new_head = node->next;
        if (new_head == nullptr) {
            // Queue is empty
            dequeue_m.unlock();
            return false;
        }
        *value = new_head->value;
        q_head_ = new_head;

        dequeue_m.unlock();
        my_allocator_.freeNode(node);
        return true;
    }

    void cleanup()
    {
        // Assume clean will only called by one thread only
        Node<T>* node = q_head_;
        while(node != nullptr) {
            Node<T>* next_node = node->next;
            my_allocator_.freeNode(node);
            node = next_node;
        }
        my_allocator_.cleanup();
    }
};
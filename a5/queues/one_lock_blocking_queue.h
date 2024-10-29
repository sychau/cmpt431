#include "../common/allocator.h"
#include <mutex>
#include <condition_variable>
#include <atomic>

template <class T>
class Node
{
public:
    T value;
    Node<T>* next;
};

// updated in both driver programs to indicate when you run out of enqueue operations
extern std::atomic<bool> no_more_enqueues; 

template <class T>
class OneLockBlockingQueue
{
    Node<T>* q_head_; // First node is sentinel
    Node<T>* q_tail_;
    std::atomic<bool> wakeup_dq; // signal for enqueuers to wake up waiting dequeuers
    CustomAllocator my_allocator_;
    std::mutex m;

public:
    OneLockBlockingQueue() : q_head_(nullptr), q_tail_(nullptr), my_allocator_() 
    {
        std::cout << "Using OneLockBlockingQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));

        // Initialize the queue head or tail here
        Node<T>* sentinel = (Node<T>*)my_allocator_.newNode();
        sentinel->value = 0;
        sentinel->next = nullptr;
        q_head_ = sentinel;
        q_tail_ = sentinel;
        wakeup_dq.store(false);
    }

    void enqueue(T value)
    {
        Node<T>* new_node = (Node<T>*)my_allocator_.newNode();
        new_node->value = value;
        new_node->next = nullptr;
        m.lock();
        q_tail_->next = new_node;
        q_tail_ = new_node;
        
        // Signal wakeup for waiting dequeue operations
        wakeup_dq.store(true);
        m.unlock();
    }

    bool dequeue(T *value)
    {
        m.lock();
        Node<T>* node = q_head_;
        Node<T>* new_head = node->next;
        int i = 0;
        while(q_head_->next == nullptr) {
            // Queue is empty
            // Wait until enqueuer wakes me up OR no more enqueues are coming
            m.unlock();
            while(!(wakeup_dq.load() || no_more_enqueues.load()));
            m.lock();
 
            if (q_head_->next != nullptr) {
                // Queue is not empty
                node = q_head_;
                new_head = q_head_->next;

                // Update wakeup_dq signal
                wakeup_dq.store(false);

            } else if (no_more_enqueues.load()) {
                // Queue empty and no more enqueues are coming
                m.unlock();
                return false;
            }
            i++;
        }

        *value = new_head->value;
        q_head_ = new_head;

        m.unlock();
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
#include "../common/allocator.h"

#define LFENCE asm volatile("lfence" : : : "memory")
#define SFENCE asm volatile("sfence" : : : "memory")

template<class P>
struct pointer_t {
    P* ptr;

    pointer_t() : ptr(nullptr) {}

    pointer_t(P* p) : ptr(p) {}

    pointer_t(uintptr_t addr, uintptr_t count) {
        ptr = (P*)((count << 48) | (addr & 0x0000FFFFFFFFFFFF));
    }

    // Credit: https://stackoverflow.com/questions/75704394
    P* address() {
        return (P*)((uintptr_t)ptr & 0x0000FFFFFFFFFFFF);
    }

    uint count() {
        return (reinterpret_cast<uintptr_t>(ptr) >> 48) & 0xFFFF;
    }

    bool operator==(const pointer_t<P>& other) const {
        return ptr == other.ptr;
    }
};

template<class T>
struct Node {
    T value;
    pointer_t<Node<T>> next;
};

template <class T>
class NonBlockingQueue
{
    pointer_t<Node<T>> q_head;
    pointer_t<Node<T>> q_tail;
    CustomAllocator my_allocator;

public:
    NonBlockingQueue() : my_allocator()
    {
        std::cout << "Using NonBlockingQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator.initialize(t_my_allocator_size, sizeof(Node<T>));

        // Initialize the queue head or tail here
        Node<T>* sentinel = (Node<T>*)my_allocator.newNode();
        sentinel->next.ptr = nullptr;

        pointer_t<Node<T>> p_sentinel(sentinel);

        q_head = p_sentinel;
        q_tail = p_sentinel;
    }

    void enqueue(T value)
    {
        // Use LFENCE and SFENCE as mentioned in pseudocode
        Node<T>* node = (Node<T>*)my_allocator.newNode();
        node->value = value;
        node->next.ptr = nullptr;
        pointer_t<Node<T>> p_node(node);
        pointer_t<Node<T>> p_tail;
        pointer_t<Node<T>> p_next;
        SFENCE;

        while (true) {
            p_tail = q_tail;
            LFENCE;
            p_next = p_tail.address()->next;
            LFENCE;
            if (p_tail == q_tail) {
                if (p_next.ptr == nullptr) {
                    if (CAS(&p_tail.address()->next, p_next,
                            pointer_t<Node<T>>(reinterpret_cast<uintptr_t>(p_node.ptr), p_next.count() + 1))) {
                        break;
                    } else {
                        CAS(&p_tail, p_tail,
                            pointer_t<Node<T>>(reinterpret_cast<uintptr_t>(p_next.ptr), p_tail.count() + 1));
                    }
                }
            }
        }
        SFENCE;
        CAS(&q_tail, p_tail, pointer_t<Node<T>>(reinterpret_cast<uintptr_t>(p_node.ptr), p_tail.count() + 1));
    }

    bool dequeue(T *value)
    {
        // Use LFENCE and SFENCE as mentioned in pseudocode
        pointer_t<Node<T>> p_head;
        pointer_t<Node<T>> p_tail;
        pointer_t<Node<T>> p_next;
        while(true){
            p_head = q_head;
            LFENCE;
            p_tail = q_tail;
            LFENCE;
            p_next = p_head.address()->next;
            LFENCE;
            if (p_head == q_head) {
                if(p_head == p_tail) {
                    if(p_next.ptr == nullptr)
                        return false;
                    CAS(&q_tail, p_tail,
                        pointer_t<Node<T>>(reinterpret_cast<uintptr_t>(p_next.ptr), p_tail.count() + 1));
                } else {
                    *value = p_next.address()->value;
                    if(CAS(&q_head, p_head,
                        pointer_t<Node<T>>(reinterpret_cast<uintptr_t>(p_next.ptr), p_head.count() + 1)))
                        break;
                }
            }
        }
        my_allocator.freeNode(p_head.address());
        return true;
    }

    void cleanup()
    {
        // Assume clean will only called by one thread only
        pointer_t<Node<T>> p_node = q_head;
        while(p_node.ptr != nullptr) {
            pointer_t<Node<T>> p_next_node = p_node.address()->next;
            my_allocator.freeNode(p_node.address());
            p_node = p_next_node;
        }
        my_allocator.cleanup();
    }
};


//
//  mymutex.h
//  searchNoise
//
//  Created by carushi on 12/04/25.
//  Copyright 2012年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_mymutex_h
#define searchNoise_mymutex_h

#include <pthread.h>

namespace mutex {
    class Mutex
    {
    private:
        pthread_mutex_t *mutex;
    public:
        Mutex(pthread_mutex_t* Mutex) : mutex(Mutex) {}
        ~Mutex() { unlock(); }
        void lock() { pthread_mutex_lock(mutex); }
        void unlock() { pthread_mutex_unlock(mutex); }
    };
}


#endif

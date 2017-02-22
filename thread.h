//
//  thread.h
//  searchNoise
//
//  Created by carushi on 12/02/16.
//  Copyright 2012年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_thread_h
#define searchNoise_thread_h

#include <iostream>
#include <vector>
#include <string>
#include "mymutex.h"

#define MAXTHREAD 8

namespace thread {

    using std::vector;
	using std::cout;
    using std::endl;
    using mutex::Mutex;

    template<class Material>
    class Thread
    {
    private:
        pthread_mutex_t *mutex;

    public:
        Thread(){}
        ~Thread(){}
        void variableThread(vector<Material> &, const int, void*(*)(void*));
    };

    template <class Material>
    void Thread<Material>::variableThread(vector<Material>& mat, const int threadnum, void*(*function)(void*))
    {
        if ( mat.size() != threadnum || threadnum > MAXTHREAD ) return;
        mutex = new pthread_mutex_t;
        pthread_mutex_init(mutex, NULL);
        pthread_t threadid[threadnum];
        vector<Mutex> mutexList = vector<Mutex> (threadnum, Mutex(mutex));
        for (int i = 0; i < threadnum; i++)
            mat[i].setMutex(&mutexList[i]);
        for (int i = 0; i < threadnum; i++)
            if ( pthread_create(&(threadid[i]), NULL, function, (void*)(&mat[i])) != 0 )
                cout << "error thread" << i << endl;
        for (int i = 0; i < threadnum; i++)
        {
            void *ret = NULL;
            if ( pthread_join(threadid[i], &ret) )
                cout << "error thread end" << i << endl;
            mat[i].setMutex(NULL);
        }
        mat.clear();
        delete mutex;
        return;
    }

}
#endif

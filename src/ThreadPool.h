
//
// ThreadPool.h
//
// Eric Viara for Institut Curie, January 2013
//
#pragma once

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <pthread.h>
#include <assert.h>
#include <iostream>

#include <vector>
#include <map>
#include <errno.h>
#include <string.h>
#include <unistd.h> //starting from v6.6 for compatibility with Ubuntu

#ifdef _WIN32
//x32 Windows definitions
#include <time.h>
#else
//other platforms
#include <sys/time.h>
#endif

class ThreadPool;
class ThreadArg;

class ThreadPoolManager {

public:
  enum Flags {
	VERBOSE = 1
  };

  static void init(unsigned int max_threads, unsigned int flags = 0) {
	assert(instance == NULL);
	instance = new ThreadPoolManager(max_threads, flags);
  }

  static ThreadPoolManager* getInstance() {
	return instance;
  }

  bool reserveOneThread();
  void releaseOneThread();
  unsigned int getMaxThreads() const {return max_threads;}
  unsigned int getFlags() const {return flags;}

  bool isMainThread() const {return pthread_self() == main_pthread;}

  ThreadPool* newThreadPool(const std::string& name);

  void lock() const;
  void unlock() const;

private:
  static ThreadPoolManager* instance;

  unsigned int max_threads;
  unsigned int flags;
  unsigned int available_threads;
  mutable pthread_mutex_t mp;
  pthread_t main_pthread;

  ThreadPoolManager(unsigned int max_threads, unsigned int flags = 0);
};

class Thread {
  static const int WAIT_USECONDS = 50000;
  typedef void* (*wrapper_t)(void*);

public:
  Thread(void* (*wrapper)(void*), ThreadArg* arg, bool delete_arg, unsigned int num) : wrapper(wrapper), arg(arg), delete_arg(delete_arg), num(num) {
	state = NOT_STARTED;
	tid = 0;
	memset(&tv_started, 0, sizeof(tv_started));
	memset(&tv_finished, 0, sizeof(tv_finished));
	pthread_mutex_init(&mp, NULL);
  }

  void running() {
	pthread_mutex_lock(&mp);
	gettimeofday(&tv_started, NULL);
	state = RUNNING;
	tid = pthread_self();
	pthread_mutex_unlock(&mp);
  }

  void finished() {
	pthread_mutex_lock(&mp);
	gettimeofday(&tv_finished, NULL);
	state = FINISHED;
	pthread_mutex_unlock(&mp);
  }

  bool waitFinished() {


   #ifdef _WIN32
//x32 Windows definitions
	Sleep(WAIT_USECONDS);
#else
//other platforms
	usleep(WAIT_USECONDS);
#endif

	usleep(WAIT_USECONDS);
	pthread_mutex_lock(&mp);
	if (state == FINISHED) {
	  pthread_mutex_unlock(&mp);
	  return true;
	}
	pthread_mutex_unlock(&mp);
	return false;
  }

  unsigned int getNumber() const {return num;}
  pthread_t getThreadSelf() const {return tid;}

  long long duration() const {
	pthread_mutex_lock(&mp);
	long long tv_usec = (state == FINISHED ? (unsigned long long)(tv_finished.tv_sec - tv_started.tv_sec) * 1000000 + (tv_finished.tv_usec - tv_started.tv_usec) : -1);
	pthread_mutex_unlock(&mp);
	return tv_usec;
  }

  ThreadArg* getArg() {return arg;}
  bool deleteArg() const {return delete_arg;}
  void* (*getWrapper())(void*) {return wrapper;}

private:
  void* (*wrapper)(void*);

  ThreadArg* arg;
  bool delete_arg;
  pthread_t tid;

  enum State {
	NOT_STARTED = 0,
	RUNNING = 1,
	FINISHED = 2
  };

  State state;
  unsigned int num;
  mutable pthread_mutex_t mp;
  struct timeval tv_started;
  struct timeval tv_finished;
};

class ThreadPool {

public:
  void addThread(void* (*wrapper)(void*), ThreadArg* arg, bool delete_arg = true) {
	thread_list.push_back(new Thread(wrapper, arg, delete_arg, thread_list.size()+1));
  }

  const std::string& getName() const {return name;}

  void run();

  ~ThreadPool();

private:
  ThreadPoolManager* thrPoolManager;
  std::string name;

  Thread* isOneThreadFinished(std::map<Thread*, bool>& thread_map);
  void checkOneThreadFinished(std::map<Thread*, bool>& thread_map);
  static void* wrapper(void* arg);

  void wait(std::map<Thread*, bool>& thread_map);

  std::vector<Thread*> thread_list;

  friend ThreadPool* ThreadPoolManager::newThreadPool(const std::string& name);

  ThreadPool(ThreadPoolManager* thrPoolManager, const std::string& name) : thrPoolManager(thrPoolManager), name(name) { }

};

class ThreadArg {

public:
  virtual ~ThreadArg() {
  }
};

#endif


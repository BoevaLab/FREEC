
//
// ThreadPool.cpp
//
// Eric Viara for Institut Curie, January 2013
//

#include <iostream>
#include "ThreadPool.h"

ThreadPoolManager* ThreadPoolManager::instance;

ThreadPoolManager::ThreadPoolManager(unsigned int _max_threads, unsigned int flags) : flags(flags)
{
  if (_max_threads == 0) {
	_max_threads = 1;
  }

  max_threads = _max_threads-1;
  available_threads = max_threads;
  pthread_mutex_init(&mp, NULL);
  main_pthread = pthread_self();
}


ThreadPool* ThreadPoolManager::newThreadPool(const std::string& name)
{
  return new ThreadPool(this, name);
}

void ThreadPoolManager::lock() const
{
  pthread_mutex_lock(&mp);
}

void ThreadPoolManager::unlock() const
{
  pthread_mutex_unlock(&mp);
}

bool ThreadPoolManager::reserveOneThread()
{
  lock();
  if (available_threads >= 1) {
	--available_threads;
	unlock();
	return true;
  }
  unlock();
  return false;
}

void ThreadPoolManager::releaseOneThread()
{
  lock();
  ++available_threads;
  unlock();
}

void ThreadPool::wait(std::map<Thread*, bool>& thread_map)
{
  while (thread_map.size() > 0) {
	checkOneThreadFinished(thread_map);
  }
}

void* ThreadPool::wrapper(void* arg)
{
  Thread* thread = (Thread*)arg;
  thread->running();
  void* status = thread->getWrapper()(thread->getArg());
  thread->finished();
  return status;
}

void ThreadPool::run()
{
  std::vector<Thread*>::iterator begin = thread_list.begin();
  std::vector<Thread*>::iterator end = thread_list.end();
  std::map<Thread*, bool> thread_map;

  bool mono_thread = thread_list.size() == 1;
  unsigned int thread_num = 1;
  time_t t0 = time(NULL);
  if (thrPoolManager->getFlags() & ThreadPoolManager::VERBOSE) {
	thrPoolManager->lock();
	std::cout << "Thread pool [" << getName() << "] started [" << thread_list.size() << " threads]\n";
	thrPoolManager->unlock();
  }

  do {
	Thread* thread = *begin++;
	pthread_t tid = 0;
	if (mono_thread) {
	  ThreadPool::wrapper(thread);
	} else {
	  pthread_create(&tid, NULL, ThreadPool::wrapper, thread);
	}
	if (thrPoolManager->getFlags() & ThreadPoolManager::VERBOSE) {
	  thrPoolManager->lock();
	  std::cout << "Thread launched [" << getName() << "#" << thread_num << ":" << tid << "] at [" << (time(NULL) - t0) << " secs]" << std::endl;
	  thrPoolManager->unlock();
	}
	thread_map[thread] = true;
	while (!thrPoolManager->reserveOneThread()) {
	  checkOneThreadFinished(thread_map);
	}
	thread_num++;
  } while (begin != end);

  wait(thread_map);

  if (thrPoolManager->getFlags() & ThreadPoolManager::VERBOSE) {
	thrPoolManager->lock();
	std::cout << "Thread pool [" << getName() << "] terminated" << std::endl;
	thrPoolManager->unlock();
  }
}

void ThreadPool::checkOneThreadFinished(std::map<Thread*, bool>& thread_map)
{
  Thread* finished_thread = isOneThreadFinished(thread_map);
  if (NULL != finished_thread) {
	if (thrPoolManager->getFlags() & ThreadPoolManager::VERBOSE) {
	  long long duration_usec = finished_thread->duration() ;
	  long long duration_sec = duration_usec / 1000000;
	  thrPoolManager->lock();
	  std::cout << "Thread terminated [" << getName() << "#" << finished_thread->getNumber() << ":" << finished_thread->getThreadSelf() << "] in [" << duration_sec << "." << (duration_usec - duration_sec*1000000) << " secs]" << std::endl;
	  thrPoolManager->unlock();
	}
	thread_map.erase(thread_map.find(finished_thread));
	thrPoolManager->releaseOneThread();
	if (finished_thread->deleteArg()) {
	  delete finished_thread->getArg();
	}
  }
}

Thread* ThreadPool::isOneThreadFinished(std::map<Thread*, bool>& thread_map)
{
  std::map<Thread*, bool>::iterator begin = thread_map.begin();
  std::map<Thread*, bool>::iterator end = thread_map.end();
  while (begin != end) {
	Thread* thread = (*begin).first;
	if (thread->waitFinished()) {
	  return thread;
	}
	++begin;
  }
  return NULL;
}

ThreadPool::~ThreadPool()
{
  std::vector<Thread*>::iterator begin = thread_list.begin();
  std::vector<Thread*>::iterator end = thread_list.end();
  while (begin != end) {
	delete *begin;
	++begin;
  }
}


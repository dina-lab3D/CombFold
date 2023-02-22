#include "Timer.h"
#include <sys/times.h>
#include <unistd.h>
#include <string.h>

Timer::Timer() : id(NULL)
{
  reset();
}

Timer::Timer(const char* const tmID)
{
  unsigned n=strlen(tmID);
  id = new char[n+1];
  strcpy(id,tmID);
  reset();
}

Timer::~Timer()
{
  if (id)
    delete [] id;
}

void Timer::reset()
{
  struct tms ptms;
  times(&ptms);

  begin_clock = ptms.tms_utime;
  clock_sum = 0;
}

void Timer::start()
{
  struct tms ptms;
  times(&ptms);

  begin_clock = ptms.tms_utime;
}

void Timer::stop()
{
  if (begin_clock >= 0) {
    struct tms ptms;
    times(&ptms);

    clock_sum+= ptms.tms_utime - begin_clock;
    begin_clock = -1;
  }
}

float Timer::operator()() const
{
  if (begin_clock >= 0){
    struct tms ptms;
    times(&ptms);

    return (clock_sum + (ptms.tms_utime - begin_clock))/(float)sysconf(_SC_CLK_TCK);//CLOCKS_PER_SEC;
  }else
    return clock_sum/(float)sysconf(_SC_CLK_TCK); //(float)CLOCKS_PER_SEC;
}

std::ostream& operator<<(std::ostream& out, const Timer& tm)
{
  if (tm.id)
    out << tm.id << ": ";
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::right, std::ios::adjustfield);
  out.precision(2);
  out.width(8);
  unsigned int time = (unsigned int)(tm()*100.0);
  unsigned int sec100 = time%100;
  time /= 100;
  unsigned int sec = time%60;
  time /= 60;
  unsigned int min = time%60;
  time /= 60;
  out.fill('0');
  out.width(2);   out << time << ':';
  out.width(2);   out << min << ':';
  out.width(2);   out << sec << '.';
  out.width(2);   out << sec100;
  out.fill(' ');
  return out;
}

#ifndef __Timer_H
#define __Timer_H

#include <ctime>
#include <string>
#include <iostream>

/*
CLASS
  Timer

  This is a general class to compute the user processor time on parts of the program.

KEYWORD
    time, clock, clocks_per_sec

AUTHORS
    Zipi Fligelman (zipo@math.tau.ac.il)

    copyright: SAMBA Group , Tel-Aviv Univ. Israel, 1998.

CHANGES LOG
<UL>
<LI>
</UL>

USAGE
    Before an important code part you define a timer thus strating it.
    After starting it every time you will output it, it will compute the
    proc sec. passed between the construction of the Timer and the curr
    time.
    You can also stop the timer after the operation you wanted has finished.
    Stopping the clock means that every time you will output the class, you
    will only recieve the time that was set between start and finish.
    You can accumulate the time on different part of the program using the
    start and stop facilities, the time will be added together. Thus giving
    you exact information on running of the program without for example major
    debugging parts

    first example:
    EXAMPLE
    Timer timer("prog_time");
    ... doing preprocessing
    ... some debugging info
    cout << timer; ( how much time preprocessing  & debugging took.)
    .... doing matching
    .... some debugging info
    cout << timer ; (how much time preprocessing , matching & debugging took.)
    .... doing verification
    cout << timer; (how much time the whole program took including debugging)
    END

    second example:
    EXAMPLE
    Timer timer("on_off_time");
    .... doing preprocessing
    timer.stop();
    .... som debugging info
    cout << timer; (how much time preprocessing took.)
    .... doing matching
    .... some debugging info
    cout << timer; (how much time preprocessing took.)
    .... doing verification
    cout << timer; (how much time preprocessing took.)
    END

    third example:
    EXAMPLE
    Timer timer("on_off_time");
    .... doing preprocessing
    timer.stop();
    .... some debugging info
    cout << timer; (how much time preprocessing took.)
    timer.start();
    .... doing matching
    cout << timer; (how much time only preprocessing and matching took.)
    timer.stop();
    .... some debugging info
    timer.start();
    .... doing verification
    cout << timer; (how much time preprocessing took.)
    END
*/
class Timer
{
public:

  // Group: Constructors

  //// Empty constructor. Resetting the start to zero. Assuming null Id.
  Timer();
  //// special possibilty if you want to keep different timers and
  // differentiate them by an ID
  Timer(const char* const tmID);

  //// destructor
  ~Timer();


  // Group: Using The class

  //// resetting the time so clock will be current clock and sum of clocks
  // will be zero
  void reset();
  //// Enable the timer to continue after a stop operation has been performed
  void start();
  //// Stopping the timer operation
  void stop();

  // Group: Operators

  //// If the timer is stopped returning the sec of proc time
  // accumalated in memory. Otherwise the time accumalated in memory and
  // the time form the last start operation until present.
  float operator()() const;

  //// printing the i.d and the time (on what time is printed see  operator())
  friend std::ostream& operator<<(std::ostream& out, const Timer& tm);

private:
  char *id;
  clock_t begin_clock, clock_sum;

};

#endif

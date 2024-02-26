#ifndef CIRCLIST_H
#define CIRCLIST_H

#include "CircularTime.h"
#include <deque>

//dir=-1 left (time increase)
//dir=1 right (time decrease)

struct Kink
{
    CircularTime ti;
    int op;
};

class CircList
{
public:
    CircList() : que() {}
    void Pop(int dir) { dir==-1 ? que.pop_front() : que.pop_back(); }
    void Push(int dir, Kink k) { dir==-1 ? que.emplace_back(k) : que.emplace_front(k); }

    Kink Top(int dir) const { return dir==-1 ? que.front() : que.back(); }
    Kink Top(int dir, int depth) const {return dir==-1 ? que[depth] : que[que.size()-depth-1];}

    std::deque<Kink>::size_type Length() {return que.size();}
    bool Empty() const {return que.empty();}
    Kink &operator[](std::deque<Kink>::size_type i) {return que[i];}

    std::deque<Kink> GetQue() { return que; }
private:
    std::deque<Kink> que;

};

#endif // CIRCLIST_H

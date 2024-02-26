#ifndef SITE_H
#define SITE_H

#include <vector>

using namespace std;

class Site
{
public:
    static int &GetSite(int x, int y)
    {
        if (x > l-1) x-=l;
        else if (x < 0) x+=l;

        if (y > l-1) y-=l;
        else if (y < 0) y+=l;

        return si[x][y];
    }

    static void SetL(int L) { l=L; }

    static void SetSite(vector<vector<int>> &s) { si=s; }

    static int L() { return l; }

    static vector<vector<int>> Si() { return si; }

private:

    static int l;

    static vector<vector<int>> si;

};

#endif // SITE_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <time.h>
#include <iomanip>
#include <bits/stdc++.h>
#include <deque>
#include <queue>
#include <unordered_map>
#include <utility>
#include <new>>
#include <string>

#define square(a) (a)*(a)

/**
      c++11
*/
using namespace std;
int N100;
int N101;
int f1, f2, f3;
bool is_reached;
int l = 0;
int bomb_path = 0;
int no_bomb_path = 0;
struct coord
{
    int x;
    int y;
    int z;
    coord(int x1, int y1, int z1)
    {
        x = x1;
        y = y1;
        z = z1;
    }
    coord()
    {
        x = -1;
        y = -1;
        z = -1;
    }
};
struct nodes_pair
{
    queue<pair<coord, bool> > pathic;
};

struct node{
    coord coordinates;
    bool is_bomb_used;
    nodes_pair path;
    //int g;
    bool is_empty()
    {
        //cout << "is_empty "<< coordinates.x << " " << coordinates.y << " " << coordinates.z << endl;
        if ((coordinates.x == -1) && (coordinates.y == -1) && (coordinates.z == -1))
            return true;
        else
            return false;
    }
    bool operator==(const node& r) const{
      return coordinates.x == r.coordinates.x && coordinates.y == r.coordinates.y && coordinates.z == r.coordinates.z;
    }

};

int failing_node;

struct node1{
    coord coordinates;
    bool is_bomb_used;
};

class MyHash{
public:
size_t operator()(const node x) const {
    return std::hash<int>()(x.coordinates.x*100+ x.coordinates.y*10+ x.coordinates.z);}
};

bool is_blackhole_zone(int x, int y, int z, bool ***blackholes)
{
    int ci[] = {-1, 0, 1};
    for (int i = 0; i<3; i++)
    {
        for(int j = 0; j<3; j++)
        {
            for (int k = 0; k<3; k++)
            {
                if((x+ci[i] >=0) && (x+ci[i]<100) && (y+ci[j]>=0) && (y+ci[j] < 100) && (z+ci[k]>=0) && (z+ci[k] < 100))
                {
                    if (blackholes[x+ci[i]][y+ci[j]][z+ci[k]])
                    {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

/*node*/void bfs_maze(bool ***black_holes, bool ***krakens, bool ***planets, int x, int y, int z)
{
 queue<node> queue1;
 int ***used;
 used = new int **[100];
 for (int i = 0; i<100; i++)
 {
     used[i] = new int *[100];
     for (int j = 0; j<100; j++)
     {
         used[i][j] = new int[100];
         for (int k = 0; k<100; k++)
         {
             used[i][j][k] = 0;
         }
     }
 }
 //used[0][0][0] = 2;
 node curr;
 curr.is_bomb_used = false;
 curr.coordinates = coord(x, y, z);
 //curr.path = NULL;
 queue1.push(curr);
 int g = 0;
 try{
 while (!queue1.empty())
 {
     g++;
     node current = queue1.front();
     queue1.pop();
     //cout << "after pop " << queue1.size() << endl;
     if(current.is_bomb_used && used[current.coordinates.x][current.coordinates.y][current.coordinates.z]!=2)
        used[current.coordinates.x][current.coordinates.y][current.coordinates.z] = 1;
     else
        used[current.coordinates.x][current.coordinates.y][current.coordinates.z] = 2;

     int xi[] = {1, 0, 0, -1, 0, 0};
     int yi[] = {0, 1, 0, 0, -1, 0};
     int zi[] = {0, 0, 1, 0, 0, -1};
     int random_seq[] = {-1, -1, -1, -1, -1, -1};
     for (int i = 0; i<6; i++)
     {
         int random_el = rand() % 5;
         //cout << "random " << random_el << endl;
         if(random_seq[random_el]!=-1)
         {
             int i = -1;
             int c = -1;
             while (i<6 && c!=-1)
             {
                 i++;
                 c = random_seq[i];
             }
             random_el = i;
             random_seq[random_el] = 1;
         }
         else
            random_seq[random_el] = 1;
         int n_x = current.coordinates.x + xi[random_el];
         int n_y = current.coordinates.y + yi[random_el];
         int n_z = current.coordinates.z + zi[random_el];
         f1 = n_x;
         f2 = n_y;
         f3 = n_z;
         if ((n_x >= 0) && (n_x < N100) && (n_y >= 0) && (n_y < N100) && (n_z >= 0) && (n_z < N100))
         {
             if((used[n_x][n_y][n_z]<2))
             {
                 if(/*!is_blackhole_zone(n_x, n_y, n_z, black_holes)*/!black_holes[n_x][n_y][n_z])
                 {
                     node tmp1;
                     tmp1 = current;
                     tmp1.coordinates.x = n_x;
                     tmp1.coordinates.y = n_y;
                     tmp1.coordinates.z = n_z;
                     //tmp1.is_bomb_used = false;
                    failing_node = tmp1.path.pathic.size();
                     if (tmp1.path.pathic.size()==0 && current.coordinates.x == 0 && current.coordinates.y == 0 && current.coordinates.z == 0)
                     {
                        tmp1.path.pathic.push(make_pair(current.coordinates, current.is_bomb_used));
                        //cout << "after this " << endl;
                     }

                     bool bomb_used = false;
                     if (planets[n_x][n_y][n_z])
                     {

                        //cout << "went through " << g << endl;
                        //tmp1.path.pathic.push(make_pair(coord(n_x, n_y, n_z), false));
                        cout << "path found by random algorithm" << endl;
                        //cout << n_x << " " << n_y << " " << n_z << endl;
                        //cout << "the size of the queue " << queue1.size() << endl;
                        cout << tmp1.path.pathic.size() << endl;
                        while (!tmp1.path.pathic.empty())
                        {
                            bool bomb = tmp1.path.pathic.front().second;
                            coord temp = tmp1.path.pathic.front().first;
                            tmp1.path.pathic.pop();
                            if (bomb)
                            {
                                cout << "M" << " " << temp.x << " " << temp.y << " " << temp.z << endl;
                            }
                            cout <<temp.x << " " << temp.y << " " << temp.z << endl;
                        }
                        //path_bomb* visited = new path_bomb;
                        //visited->is_bomb_used = false;
                        cout << tmp1.coordinates.x << " " << tmp1.coordinates.y << " " << tmp1.coordinates.z << endl;
                        return;
                        // return tmp1;
                     }
                     if (krakens[n_x][n_y][n_z])
                     {
                         if(current.is_bomb_used)
                         {

                             continue;
                         }
                         tmp1.is_bomb_used= true;
                         bomb_used = true;
                     }
                     tmp1.path.pathic.push(make_pair(tmp1.coordinates, bomb_used));
                     /**********************************
                     ****TO REMOVE**********************/
                      if(tmp1.is_bomb_used && used[tmp1.coordinates.x][tmp1.coordinates.y][tmp1.coordinates.z]!=2)
                            used[tmp1.coordinates.x][tmp1.coordinates.y][tmp1.coordinates.z] = 1;
                         else
                            used[tmp1.coordinates.x][tmp1.coordinates.y][tmp1.coordinates.z] = 2;
                    /***********************************
                    ***********************************/
                     //cout << n_x << " " << n_y << " " << n_z << " ";
                     queue1.push(tmp1);
                 }
             }
         }
     }
     //cout << "queue1 after push " << queue1.size() << endl;
     //cout << " \n" << endl;

 }
 }
 catch(const std::bad_alloc& e)
 {
     std::cout << "It took to much stack space to store the values in the queue. Allocation failed when adding the node: " << f1 << " - " << " " << f2 << " - " << " " << f3 << " - " << " queue1 size " << queue1.size() << e.what() << '\n';
     cout << "size of the tmp1 node " << failing_node << endl;
 }
 cout << "This is unsolvable mission." << endl;
}
struct path_bomb
{
    deque<node1> path;
    bool is_bomb_used;
};
int xi1[] = {1, 0, 0, -1, 0, 0};
int yi1[] = {0, 1, 0, 0, -1, 0};
int zi1[] = {0, 0, 1, 0, 0, -1};
path_bomb* found2 = new path_bomb;
bool backtracking(bool ***black_holes, bool ***kracens, bool ***planets, int ***used, int x, int y, int z, path_bomb  *visited, bool is_bomb_used, int p_count, int counts = 0)
{
    if (p_count == 0)
        return false;
    node1 new_v;
    new_v.coordinates = coord(x, y, z);
    new_v.is_bomb_used = false;
    bool flag = false;

    if (planets[x][y][z])
    {
        visited->path.push_back(new_v);
        int j = 0;
        cout << visited->path.size() << endl;
        std::cin.ignore();
        while(!visited->path.empty())
        {
            j++;
            if (visited->path.front().is_bomb_used)
        {
            cout << "M" << " " << visited->path.front().coordinates.x << " " << visited->path.front().coordinates.y << " " << visited->path.front().coordinates.z << endl;
        }
            cout << visited->path.front().coordinates.x << " " <<
            visited->path.front().coordinates.y << " " <<
            visited->path.front().coordinates.z << endl;
            visited->path.pop_front();
        }
        cout << "the number of nodes is " << j << endl;
        //std::cin.ignore();
        return true;
    }
        if(/*!is_blackhole_zone(x, y, z, black_holes)*/!black_holes[x][y][z])
        {
             if (kracens[x][y][z] && !visited->is_bomb_used)
             {
                 visited->is_bomb_used = true;
                 flag = true;
              }
             else if (!kracens[x][y][z])
             {
                flag = true;
             }

         if (flag)
         {
        /*else see if it is safe to make a move*/
            /*********/
            cout << "here" << endl;
            cout << x << " " << y << " " << z << endl;
            cout<< visited->is_bomb_used << endl;
            visited->path.push_back(new_v);

            for (int i = 0; i<6; i++)
            {
                int n_x = x + xi1[i];
                int n_y = y + yi1[i];
                int n_z = z + zi1[i];
                cout << n_x << " " << n_y << " " << n_z << endl;
                if ((n_x >= 0) && (n_x < 100) && (n_y >= 0) && (n_y < 100) && (n_z >= 0) && (n_z < 100
                                                                                             ))
                {
                    if(used[n_x][n_y][n_z]<2)
                    {
                        if(visited->is_bomb_used && used[n_x][n_y][n_z]==0)
                        {
                            used[n_x][n_y][n_z] = 1;
                        }
                       else if(visited->is_bomb_used && (used[n_x][n_y][n_z]==1 || used[n_x][n_y][n_z]==0))
                            used[n_x][n_y][n_z] = 2;
                        else continue;
                        if (!black_holes[n_x][n_y][n_z])
                        {
                            if (kracens[n_x][n_y][n_z] && !visited->is_bomb_used)
                            {
                                 cout << n_x << " " << n_y << " " << n_z << endl;
                                if(backtracking(black_holes, kracens, planets, used, n_x, n_y, n_z, visited, is_bomb_used, p_count, counts+1))
                                    return true;
                            }
                            else if(!kracens[n_x][n_y][n_z])
                            {
                                 cout << n_x << " " << n_y << " " << n_z << endl;
                                if(backtracking(black_holes, kracens, planets, used, n_x, n_y, n_z, visited, is_bomb_used, p_count, counts+1))
                                    return true;
                            }
                        }
                    }

                }
            }
            node1 del = visited->path.back();
            if (kracens[del.coordinates.x][del.coordinates.y][del.coordinates.z])
            {
                visited->is_bomb_used = false;
            }
            visited->path.pop_back();
            }
    }

    return false;
}

struct PriorityQueue
{
    typedef pair<double, node> element;
    class Prioritize{
        public:
            int operator()(const element& p1,
                           const element& p2)
                           {
                               return p1.first > p2.first;
                           }
    };
    std::priority_queue<element, std::vector<element>,
    Prioritize> elements;
    inline bool empty() const
    {
        return elements.empty();
    }

    inline void put(node item, double priority)
    {
        element el1 = std::make_pair(priority, item);
        elements.emplace(el1);
    }

    inline node get()
    {
        node best = elements.top().second;
        elements.pop();
        return best;
    }
};

bool valid(int x, int y, int z)
{
    if (x>=0 && x<N100 && y>=0 && y<N100 && z>=0 && z<N100)
    {
        return true;
    }
    else
        return false;
}

std::list<node> planets_list;

inline double heuristic(node a)
{

    int x1, y1, z1;
    x1 = a.coordinates.x;
    y1 = a.coordinates.y;
    z1 = a.coordinates.z;
    double min1 = abs(x1 - planets_list.front().coordinates.x) + abs(y1 - planets_list.front().coordinates.y) + abs(z1 - planets_list.front().coordinates.z);
    node goal2;
    for(std::list<node>::iterator it = planets_list.begin(); it != planets_list.end(); it++)
    {
    int k = abs(x1 - (*it).coordinates.x) + abs(y1 - (*it).coordinates.y) + abs(z1 - (*it).coordinates.z);
         if (k<min1)
            min1 = k;
            goal2 = (*it);
    }

    int dx1 = x1 - goal2.coordinates.x;
    int dy1 = y1 - goal2.coordinates.y;
    int dz1 = z1 - goal2.coordinates.z;
    int dx2 = 0 - goal2.coordinates.x;
    int dy2 = 0 - goal2.coordinates.y;
    int dz2 = 0 - goal2.coordinates.z;
    double cross = sqrt(square(dy1*dz2 - dz1*dy2) + square(dz1*dx2 - dx1*dz2) + square(dx1*dy2 - dy1*dx2));
    min1 = min1 + cross*0.001;
    return min1;
}

node goal;
int size_a_star = 0;
void reconstruct_path(node start, node goal, std::unordered_map<node, pair<node, node>, MyHash>*& came_from)
{
    vector<node> path;
    node current = goal;
    path.push_back(current);
    int k = 0;
    cout << bomb_path << endl;
    cout << no_bomb_path << endl;
    if ((bomb_path<no_bomb_path && bomb_path != 0)|| no_bomb_path==0)
    {
        while(!current.is_bomb_used)
        {
         current = came_from->at(current).second;
         path.push_back(current);
         size_a_star++;
        }
        current = came_from->at(current).second;
        path.push_back(current);
        size_a_star++;
        while(!(current == start))
        {
            current = came_from->at(current).first;
            path.push_back(current);
            size_a_star++;
        }
    }
    else
    {
        while(!(current == start))
        {
            current = came_from->at(current).first;
            path.push_back(current);
            size_a_star++;
        }
    }

    std::reverse(path.begin(), path.end());
    cout << "the path found by the a*" << endl;
    if (size_a_star==0)
    {
        cout << "This is unsolvable mission" << endl;
    }
    else
        {
            cout<<size_a_star << endl;
            for(node n: path)
            {
                if (n.is_bomb_used)
                {
                    cout << "M" << " " << n.coordinates.x << " " << n.coordinates.y << " " << n.coordinates.z << endl;
                }
                cout << n.coordinates.x << " " << n.coordinates.y << " " << n.coordinates.z << endl;
            }
        }
        //return path;
}

void a_star(bool ***black_holes, bool ***krakens, bool ***planets, int*** used, int x, int y, int z, unordered_map<node, pair<node, node>, MyHash> *&came_from, unordered_map<node, double, MyHash> cost, unordered_map<node, double, MyHash> bomb_cost, int p)
{
    PriorityQueue frontier;
    node start;
    start.coordinates = coord(x, y, z);
    start.is_bomb_used = false;
    bool really_bomb_used = false;
    frontier.put(start, 0);
    node no_data;
    came_from->insert(std::make_pair(start, make_pair(start, no_data)));
    cost[start] = 0;
    used[0][0][0] = 2;
    if (p==0)
        return;

    while (!frontier.empty())
    {
        node current = frontier.get();
        if((!came_from->at(current).second.is_empty()) && used[current.coordinates.x][current.coordinates.y][current.coordinates.z]!=2)
            used[current.coordinates.x][current.coordinates.y][current.coordinates.z] = 1;
        else
            used[current.coordinates.x][current.coordinates.y][current.coordinates.z] = 2;
        int xi[] = {1, 0, 0, -1, 0, 0};
        int yi[] = {0, 1, 0, 0, -1, 0};
        int zi[] = {0, 0, 1, 0, 0, -1};
        for (int i = 0; i<6; i++)
        {
            int n_x = current.coordinates.x + xi[i];
            int n_y = current.coordinates.y + yi[i];
            int n_z = current.coordinates.z + zi[i];
            if (valid(n_x, n_y, n_z))
            {
                if((used[n_x][n_y][n_z]<2))
                {
                        /** checks **/
                        node next;
                        next.coordinates = coord(n_x, n_y, n_z);
                        next.is_bomb_used = false;
                        double new_cost_bomb;
                        double new_cost;
                        if (bomb_cost.count(current)&&(!came_from->at(current).second.is_empty()))
                            new_cost_bomb = bomb_cost[current] + 1;
                        if (!came_from->at(current).first.is_empty())
                            new_cost = cost[current] + 1;
                        bool flag = false;
                        if (planets[n_x][n_y][n_z])
                            {
                                goal.coordinates = coord(n_x, n_y, n_z);
                                if (!came_from->at(current).second.is_empty())
                                {
                                    bomb_path = bomb_cost[current] + 1;
                                    came_from->insert({next, std::make_pair(no_data, current)});
                                }
                                if (!came_from->at(current).first.is_empty())
                                {
                                    if(came_from->count(next)&&!came_from->at(next).second.is_empty())
                                    {
                                        came_from->at(next) = std::make_pair(current, came_from->at(next).second);
                                    }
                                    else
                                    {
                                         came_from->insert({next, std::make_pair(current, no_data)});
                                    }
                                    no_bomb_path = cost[current] + 1;
                                }
                                is_reached=true;
                                return;
                            }

                        if(!black_holes[n_x][n_y][n_z])
                            {
                                if (krakens[n_x][n_y][n_z])
                                     {
                                         if (came_from->at(current).first.is_empty()&&(!came_from->at(current).second.is_empty()))
                                         {
                                             continue;
                                         }
                                         flag = true;
                                         new_cost_bomb = cost[current] + 1;
                                         next.is_bomb_used = true;
                                     }
                                        if(((!bomb_cost.count(next)||new_cost_bomb <bomb_cost[next]) || (!cost.count(next) || new_cost < cost[next])))
                                        {
                                            double priority;
                                            if (!came_from->count(next))
                                            {

                                                if(!flag && !came_from->at(current).first.is_empty() && !came_from->at(current).second.is_empty() && (!cost.count(next) || (new_cost<cost[next] || cost[next] == 0)) && (bomb_cost.count(next) || bomb_cost[next] == 0 || bomb_cost[next]>=new_cost_bomb))
                                                {
                                                    cost[next] = new_cost;
                                                    bomb_cost[next] = new_cost_bomb;
                                                    if(new_cost<new_cost_bomb)
                                                        priority = new_cost + heuristic(next);
                                                    else
                                                        priority = new_cost_bomb + heuristic(next);

                                                    frontier.put(next, priority);
                                                    came_from->insert({next, make_pair(current, current)});
                                                }
                                                else
                                                if (!flag && !came_from->at(current).first.is_empty() && (!cost.count(next) || new_cost<cost[next] || cost[next] == 0))
                                                {
                                                    cost[next] = new_cost;
                                                    priority = new_cost + heuristic(next);
                                                    frontier.put(next, priority);
                                                    came_from->insert({next, make_pair(current, no_data)});
                                                }
                                                else
                                                if (!flag && !came_from->at(current).second.is_empty() && (!bomb_cost.count(next) || bomb_cost[next]==0 || new_cost_bomb<bomb_cost[next]))
                                                {
                                                    bomb_cost[next] = new_cost_bomb;
                                                    priority =new_cost_bomb + heuristic(next);
                                                    frontier.put(next, priority);
                                                     came_from->insert({next, make_pair(no_data, current)});
                                                }
                                                    else if(flag && !came_from->at(current).first.is_empty() && (!bomb_cost.count(next) || bomb_cost[next] == 0 || bomb_cost[next]>new_cost_bomb))
                                                        {
                                                                bomb_cost[next] = new_cost_bomb;
                                                                priority = new_cost_bomb + heuristic(next);
                                                                frontier.put(next, priority);
                                                                came_from->insert({next, make_pair(no_data, current)});
                                                        }
                                                        else
                                                            continue;

                                            }
                                            else
                                            {
                                              if (!flag && !came_from->at(current).second.is_empty() && !came_from->at(current).first.is_empty() && (!cost.count(next) || cost[next] == 0 || new_cost<cost[next]) && (!bomb_cost.count(next) ||  bomb_cost[next] ==0 || bomb_cost[next]>new_cost_bomb))
                                              {
                                                    cost[next] = new_cost;
                                                    bomb_cost[next] = new_cost_bomb;
                                                    if(new_cost<new_cost_bomb)
                                                        priority = new_cost + heuristic(next);
                                                    else
                                                        priority = new_cost_bomb + heuristic(next);

                                                    frontier.put(next, priority);
                                                    came_from->at(next) = make_pair(current, current);
                                              }

                                                else if (!flag && !came_from->at(current).second.is_empty() && (!bomb_cost.count(next)|| bomb_cost[next] == 0 || bomb_cost[next]>new_cost_bomb))
                                              {
                                                    cost[next] = new_cost;
                                                    priority = new_cost + heuristic(next);
                                                    frontier.put(next, priority);
                                                    came_from->at(next)= make_pair(came_from->at(next).first, current);
                                              }
                                              else if (!flag && !came_from->at(current).first.is_empty() && (!cost.count(next) || cost[next] == 0 || cost[next]>new_cost))
                                              {
                                                    bomb_cost[next] = new_cost_bomb;
                                                    priority =new_cost_bomb + heuristic(next);
                                                    frontier.put(next, priority);
                                                    came_from->at(next)= make_pair(current, came_from->at(next).second);
                                              }
                                               else
                                              if( flag&&(!came_from->at(current).first.is_empty()) && (!bomb_cost.count(next) || bomb_cost[next] == 0 || bomb_cost[next]>new_cost_bomb))
                                              {
                                                    bomb_cost[next] = new_cost_bomb;
                                                    priority = new_cost_bomb + heuristic(next);
                                                    frontier.put(next, priority);
                                                    came_from->at(next)= make_pair(came_from->at(next).first, current);
                                              }
                                            }
                                            if(!came_from->at(next).second.is_empty() && used[next.coordinates.x][next.coordinates.y][next.coordinates.z]!=2)
                                                used[next.coordinates.x][next.coordinates.y][next.coordinates.z] = 1;
                                             else
                                                used[next.coordinates.x][next.coordinates.y][next.coordinates.z] = 2;
                                        }
                            }
                }
            }
        }
    }
    int k = came_from->size();
    cout << "The mission is unsolvable : -/ queue size " << frontier.elements.size() << " number of nodes in the map " << k << " the number of adds " << l << endl;
}

int main(int )
{
    N101 = 101;
    N100 = 100;
    char Type;
    int x, y, z;
    std::ofstream myfile;
    node n;
    bool ***black_holes;
    black_holes = new bool **[N101];
    bool ***krakens;
    krakens = new bool **[N101];
    bool ***planets;
    planets = new bool **[N101];
    bool ***used;
    used = new bool **[N101];
    int ***used1;
    used1 = new int **[N101];
    for (int i = 0; i<=N100; i++)
        {
            used1[i] = new int *[N101];
            for (int j = 0; j<=N100; j++)
            {
                used1[i][j] = new int[N101];
                for (int k = 0; k<=N100; k++)
                {
                    used1[i][j][k] = 0;
                }
            }
        }
    /***
    RUN TESTS ONE BY ONE
    don't use the loop like for 0 to 16, cause the names of the test files are different
    there is no test input(11).txt for example
    ***/
    for (int i = 4; i<=4; i++)
    {
            stringstream convert;
            convert << i;
            string k = convert.str();
           // myfile.open("input("+k+").txt");
            std::ifstream infile("input("+k+").txt");
            cout << "       The test " << i << " is running\n" << endl;

            int bh_count = 0;
            int k_count = 0;
            int p_count = 0;
            /**
            >>>>>>>> optimaze
            */
            for (int i = 0; i<=N100; i++)
            {
                black_holes[i] = new bool *[N101];
                krakens[i] = new bool *[N101];
                planets[i] = new bool *[N101];
                used[i] = new bool *[N101];
                for (int j = 0; j<=N100; j++)
                {
                    black_holes[i][j] = new bool [N101];
                    krakens[i][j] = new bool [N101];
                    planets[i][j] = new bool [N101];
                    used[i][j] = new bool[N101];
                    for (int k = 0; k<=N100; k++)
                    {
                        black_holes[i][j][k] = false;
                        krakens[i][j][k] = false;
                        planets[i][j][k] = false;
                        used[i][j][k] = false;
                    }
                }
            }

            while(infile >> Type >> x >> y >> z)
            {
                if (Type == 'B')
                {
                    black_holes[x][y][z] = true;
                    bh_count++;
                }
                else
                    if (Type == 'K')
                {
                    krakens[x][y][z] = true;
                    k_count++;
                }
                else
                    if (Type == 'P')
                {
                    node node_p;
                    node_p.coordinates = coord(x, y, z);
                    planets_list.push_back(node_p);
                    planets[x][y][z] = true;
                    p_count++;
                }
            }

            clock_t end;
            clock_t begin;
            is_reached = false;
            bool ***planets2;
            bool ***krakens2;
            bool ***black_holes2;
            //used
            /**
            *************RANDOM***********
            */
            begin = clock();
            /*node found = */bfs_maze(black_holes, krakens, planets, 0, 0, 0);
            end = clock();

            double print = (end - begin);
            print = print;
            printf("%.4f msec\n", print);
            cout << "end" << endl;


            deque<node1>* b = new deque<node1>;


            found2->is_bomb_used = false;

            used1[0][0][0] = 2;
            int ***used2;
            used2 = used1;
            cout << "\n" << endl;
            /************************************
            ************ A-Star *****************
            *************************************
            For A* algorithm it will be efficient
            to use priority queue
            ************************************/
            unordered_map<node, double, MyHash> cost;
            unordered_map<node, double, MyHash> bomb_cost;
            unordered_map<node, pair<node, node>, MyHash> *came_from = new unordered_map<node, pair<node, node>, MyHash>;
            begin = clock();
            a_star(black_holes, krakens, planets, used2, 0, 0, 0, came_from, cost, bomb_cost, p_count);
            end = clock();
            node start;
            start.coordinates = coord(0, 0, 0);
            start.is_bomb_used = false;
            if(is_reached)
            {
            reconstruct_path(start, goal, came_from);
            }
            print = (end - begin);
            print = print;
            printf("%.4f msec\n", print);
            cout << "\n" << endl;
            /************************
    *************************
    ******BACKTRACKING*******
    *************************
    #########################
    ****TO RUN UNCOMMENT*****
    #########################
    !!!!!!!!!!!!!!!!!!!!!!!!!
    run it carefully, needed
    to increase the stack size
    1) project->build options->
    ->linker settings->other linker options
    2) -Wl,--stack,52800000
    3) Build->Rebuild
    !!!!!!!!!!!!!!!!!!!!!!!!!
    #########################
    ***************************/

   path_bomb *visited = new path_bomb();
    visited->is_bomb_used = false;
            for (int i = 0; i<=N100; i++)
            {
                for (int j = 0; j<=N100; j++)
                {
                    for (int k = 0; k<=N100; k++)
                    {
                        used1[i][j][k] = 0;
                    }
                }
            }
            used1[0][0][0] = 2;
    cout << " \n" << endl;
    cout << "path found by backtracking" << endl;
    begin = clock();
    bool t = backtracking(black_holes, krakens, planets, used1, 0, 0, 0, visited, false, p_count);
    end = clock();
    if (!t)
        cout << "This is unsolvable mission" << endl;
    print = (end - begin);
    print = print;
    printf("%.4f msec\n", print);
    cout << "\n" << endl;
    }
    cout << "Hello world!" << endl;
    return 0;
}



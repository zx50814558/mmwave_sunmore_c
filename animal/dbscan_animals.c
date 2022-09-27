/* Copyright 2015 Gagarine Yaikhom (MIT License) */
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define UNCLASSIFIED -1
#define NOISE -2

#define CORE_POINT 1
#define NOT_CORE_POINT 0

#define SUCCESS 0
#define FAILURE -3

typedef struct point_s point_t;
struct point_s {
    double x, y, z;
    int cluster_id;
};

typedef struct node_s node_t;
struct node_s {
    unsigned int index;
    node_t *next;
};

typedef struct epsilon_neighbours_s epsilon_neighbours_t;
struct epsilon_neighbours_s {
    unsigned int num_members;
    node_t *head, *tail;
};

node_t *create_node(unsigned int index);
int append_at_end(
     unsigned int index,
     epsilon_neighbours_t *en);
epsilon_neighbours_t *get_epsilon_neighbours(
    unsigned int index,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    double (*dist)(point_t *a, point_t *b));

void destroy_epsilon_neighbours(epsilon_neighbours_t *en);
void dbscan(
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b));
int expand(
    unsigned int index,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b));
int spread(
    unsigned int index,
    epsilon_neighbours_t *seeds,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b));
double euclidean_dist(point_t *a, point_t *b);
double adjacent_intensity_dist(point_t *a, point_t *b);
unsigned int parse_input(
    FILE *file,
    point_t **points,
    double *epsilon,
    unsigned int *minpts);
unsigned int parse_input_1(
    point_t **points,
    float (*pos1X)[6],
    double *epsilon,
    unsigned int *minpts,
    int number_of_points);
void print_points(
    point_t *points,
    unsigned int num_points);

node_t *create_node(unsigned int index)
{
    node_t *n = (node_t *) calloc(1, sizeof(node_t));
    if (n == NULL)
        perror("Failed to allocate node.");
    else {
        n->index = index;
        n->next = NULL;
    }
    return n;
}

int append_at_end(
     unsigned int index,
     epsilon_neighbours_t *en)
{
    node_t *n = create_node(index);
    if (n == NULL) {
        free(en);
        return FAILURE;
    }
    if (en->head == NULL) {
        en->head = n;
        en->tail = n;
    } else {
        en->tail->next = n;
        en->tail = n;
    }
    ++(en->num_members);
    return SUCCESS;
}

epsilon_neighbours_t *get_epsilon_neighbours(
    unsigned int index,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    double (*dist)(point_t *a, point_t *b))
{
    epsilon_neighbours_t *en = (epsilon_neighbours_t *)
        calloc(1, sizeof(epsilon_neighbours_t));
    if (en == NULL) {
        perror("Failed to allocate epsilon neighbours.");
        return en;
    }
    for (int i = 0; i < num_points; ++i) {
        if (i == index)
            continue;
        if (dist(&points[index], &points[i]) > epsilon)
            continue;
        else {
            if (append_at_end(i, en) == FAILURE) {
                destroy_epsilon_neighbours(en);
                en = NULL;
                break;
            }
        }
    }
    return en;
}

void destroy_epsilon_neighbours(epsilon_neighbours_t *en)
{
    if (en) {
        node_t *t, *h = en->head;
        while (h) {
            t = h->next;
            free(h);
            h = t;
        }
        free(en);
    }
}

void dbscan(
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b))
{
    unsigned int i, cluster_id = 0;
    for (i = 0; i < num_points; ++i) {
        if (points[i].cluster_id == UNCLASSIFIED) {
            if (expand(i, cluster_id, points,
                       num_points, epsilon, minpts,
                       dist) == CORE_POINT)
                ++cluster_id;
        }
    }
}

int expand(
    unsigned int index,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b))
{
    int return_value = NOT_CORE_POINT;
    epsilon_neighbours_t *seeds =
        get_epsilon_neighbours(index, points,
                               num_points, epsilon,
                               dist);
    if (seeds == NULL)
        return FAILURE;

    if (seeds->num_members < minpts)
        points[index].cluster_id = NOISE;
    else {
        points[index].cluster_id = cluster_id;
        node_t *h = seeds->head;
        while (h) {
            points[h->index].cluster_id = cluster_id;
            h = h->next;
        }

        h = seeds->head;
        while (h) {
            spread(h->index, seeds, cluster_id, points,
                   num_points, epsilon, minpts, dist);
            h = h->next;
        }

        return_value = CORE_POINT;
    }
    destroy_epsilon_neighbours(seeds);
    return return_value;
}

int spread(
    unsigned int index,
    epsilon_neighbours_t *seeds,
    unsigned int cluster_id,
    point_t *points,
    unsigned int num_points,
    double epsilon,
    unsigned int minpts,
    double (*dist)(point_t *a, point_t *b))
{
    epsilon_neighbours_t *spread =
        get_epsilon_neighbours(index, points,
                       num_points, epsilon,
                       dist);
    if (spread == NULL)
        return FAILURE;
    if (spread->num_members >= minpts) {
        node_t *n = spread->head;
        point_t *d;
        while (n) {
            d = &points[n->index];
            if (d->cluster_id == NOISE ||
                d->cluster_id == UNCLASSIFIED) {
                if (d->cluster_id == UNCLASSIFIED) {
                    if (append_at_end(n->index, seeds)
                        == FAILURE) {
                        destroy_epsilon_neighbours(spread);
                        return FAILURE;
                    }
                }
                d->cluster_id = cluster_id;
            }
            n = n->next;
        }
    }

    destroy_epsilon_neighbours(spread);
    return SUCCESS;
}

double euclidean_dist(point_t *a, point_t *b)
{
    return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2) + pow(a->z - b->z, 2));
}
/*
void input_data(float ** intput_array)
{

}
*/
unsigned int parse_input(FILE *file, point_t **points, double *epsilon, unsigned int *minpts)
{
    unsigned int num_points, i = 0;

    //fscanf(file, "%lf %u %u\n", epsilon, minpts, &num_points);
    *epsilon = 1;
    *minpts = 2;
    num_points = 3;
    //解決第一行輸入問題 epsilon minpts num_points
    //printf("epsilon =  %lf \n", *epsilon);
    //printf("minpts =  %d \n", *minpts);
    //printf("num_points =  %d \n", num_points);
    point_t *p = (point_t *)
    //printf("*pointer 的值：%d\n", *pointer);
    calloc(num_points, sizeof(point_t)); //動態分配記憶體
    if (p == NULL) {
        perror("Failed to allocate points array");
        return 0;
    }
    /*
    while (i < num_points) {
          fscanf(file, "%lf %lf %lf\n",
                 &(p[i].x), &(p[i].y), &(p[i].z));
          p[i].cluster_id = UNCLASSIFIED;
          ++i;
    }
    */
    //load
    float pos1X[3][4];
    pos1X[0][0] = 1.0;
    pos1X[0][1] = 15.0;
    pos1X[0][2] = 15.0;
    pos1X[0][3] = -1.0;
    pos1X[1][0] = 1.0;
    pos1X[1][1] = 1.0;
    pos1X[1][2] = 1.0;
    pos1X[1][3] = -1.0;
    pos1X[2][0] = 2.0;
    pos1X[2][1] = 2.0;
    pos1X[2][2] = 12.0;
    pos1X[2][3] = -1.0;
    //只要傳入pos1X
    int row = sizeof(pos1X) / sizeof(pos1X[0]);
    printf("Number of rows: %d\n", row);
    for(int num=0; num < row; ++num)
    {
        p[num].x = pos1X[num][0];
        p[num].y = pos1X[num][1];
        p[num].z = pos1X[num][2];
        p[num].cluster_id = pos1X[num][3];
    }
    *points = p;
    return num_points;
}

void print_points(point_t *points, unsigned int num_points)
{
    unsigned int i = 0;
    
    printf("Number of points: %u\n"
        " x     y     z     cluster_id\n"
        "-----------------------------\n"
        , num_points);
        
    while (i < num_points) {
          printf("%5.2lf %5.2lf %5.2lf: %d\n", points[i].x, points[i].y, points[i].z, points[i].cluster_id);
          ++i;
    }
}

struct set_point_clude {
    float vsos[1000][4];
};
typedef struct set_point_clude Struct;
  
Struct dbscan_output(float (*pos1X)[6], int number_of_points)
{
    Struct output;
    
    point_t *points;
    double epsilon;
    unsigned int minpts;
    //unsigned int num_points =parse_input_1(stdin, &points, &epsilon, &minpts);
    unsigned int num_points =parse_input_1(&points, pos1X, &epsilon, &minpts, number_of_points);
    //point都填入了
    if (num_points) {
        dbscan(points, num_points, epsilon,minpts, euclidean_dist);
        print_points(points, num_points);
    }

    int num_points_1;
    //num_points_1 = (int ) num_points;
    num_points_1 = number_of_points;
    float dbscan_array[num_points_1][4];
    //printf("number_of_points = %d\n", number_of_points);
    //printf("%d\n", points[0].cluster_id);
    for(int num=0; num < num_points_1; ++num)
    {
        output.vsos[num][0]=points[num].x;
        output.vsos[num][1]=points[num].y;
        output.vsos[num][2]=points[num].z;
        output.vsos[num][3]=(float) points[num].cluster_id;
    }
    free(points);
    return output;//增加的
}
unsigned int parse_input_1(point_t **points, float (*pos1X)[6], double *epsilon, unsigned int *minpts , int number_of_points )
{
    unsigned int num_points, i = 0;
    *epsilon = 0.5;
    *minpts = 18;
    num_points = number_of_points;
    point_t *p = (point_t *)
    calloc(num_points, sizeof(point_t)); //動態分配記憶體
    //只要傳入pos1X
    //printf("pos1X[0][0]=%f\n", pos1X[0][0]);
    int row = number_of_points;
    for(int num=0; num < row; ++num)
    {
        p[num].x = pos1X[num][0];
        p[num].y = pos1X[num][1];
        p[num].z = pos1X[num][2];
        p[num].cluster_id = pos1X[num][3];
    }
    *points = p;
    return num_points;
    
}

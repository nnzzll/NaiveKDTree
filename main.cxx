#include "KDTree.h"
#include <algorithm>
#include <array>
#include <vector>
#include <iostream>
#include <chrono>

#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPointSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkNew.h>

template <typename T>
struct Point3D
{
    T x, y, z;
    int index;
    Point3D() : x(0), y(0), z(0), index(-1){};
    Point3D(T a, T b, T c) : x(a), y(b), z(c), index(-1){};
    Point3D(T a, T b, T c, int idx) : x(a), y(b), z(c), index(idx){};
    inline T &operator[](int i) { return i == 0 ? x : i == 1 ? y
                                                             : z; };
};

template <typename T>
struct Point2D
{
    T x, y;
    int index;
    Point2D() : x(0), y(0), index(-1){};
    Point2D(T a, T b) : x(a), y(b), index(-1){};
    Point2D(T a, T b, int idx) : x(a), y(b), index(idx){};

    inline T &operator[](int i) { return i == 0 ? x : y; };
};

int main()
{
    int N = 500;
    // Create some random points
    vtkNew<vtkPointSource> pointSource;
    pointSource->SetNumberOfPoints(N);
    pointSource->Update();

    std::vector<Point3D<double>> datasets;
    vtkPoints *randPts = pointSource->GetOutput()->GetPoints();
    for (vtkIdType i = 0; i < N; i++)
    {
        double pts[3];
        randPts->GetPoint(i, pts);
        datasets.push_back(Point3D<double>(pts[0], pts[1], pts[2], i));
        // std::cout << pts[0] << "," << pts[1] << "," << pts[2] << std::endl;
    }

    auto t1 = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::system_clock::now().time_since_epoch())
                  .count();
    
    // Create the tree
    vtkNew<vtkKdTreePointLocator> pointTree;
    pointTree->SetDataSet(pointSource->GetOutput());
    pointTree->BuildLocator();

    // Find the k closest points to (0,0,0)
    unsigned int k = 1;
    vtkNew<vtkPointSource> testSource;
    testSource->SetNumberOfPoints(1);
    testSource->Update();
    double testPoint[3];
    testSource->GetOutput()->GetPoints()->GetPoint(0, testPoint);
    vtkNew<vtkIdList> result;
    std::cout << "Test Point: " << testPoint[0] << "," << testPoint[1] << "," << testPoint[2] << std::endl;

    pointTree->FindClosestNPoints(k, testPoint, result);

    for (vtkIdType i = 0; i < k; i++)
    {
        vtkIdType point_ind = result->GetId(i);
        double p[3];
        pointSource->GetOutput()->GetPoint(point_ind, p);
        std::cout << "Closest point " << i << ": Point " << point_ind << ": ("
                  << p[0] << ", " << p[1] << ", " << p[2] << ")" << std::endl;
    }
    auto t2 = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::system_clock::now().time_since_epoch())
                  .count();

    // Should return:
    // Closest point 0: Point 2: (-0.136162, -0.0276359, 0.0369441)

    // std::vector<Point2D<double>> datasets = {Point2D<double>(7, 2, 0),
    //                                       Point2D<double>(5, 4, 1),
    //                                       Point2D<double>(9, 6, 2),
    //                                       Point2D<double>(2, 3, 3),
    //                                       Point2D<double>(4, 7, 4),
    //                                       Point2D<double>(8, 1, 5)};
    KDTree<Point3D<double>> tree(datasets, 3);
    // tree.Print();
    std::cout << tree.findNearestPoint(Point3D<double>(testPoint[0], testPoint[1], testPoint[2])) << std::endl;
    auto t3 = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::system_clock::now().time_since_epoch())
                  .count();
    std::cout << "VTK Time:" << t2 - t1 << " ms" << std::endl;
    std::cout << "MY Time:" << t3 - t2 << " ms" << std::endl;
    return EXIT_SUCCESS;
}
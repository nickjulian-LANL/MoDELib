/*****************************************************************************/
/*                                                                           */
/*  (tricall.c)                                                              */
/*                                                                           */
/*  Example program that demonstrates how to call Triangle.                  */
/*                                                                           */
/*  Accompanies Triangle Version 1.6                                         */
/*  July 19, 1996                                                            */
/*                                                                           */
/*  This file is placed in the public domain (but the file that it calls     */
/*  is still copyrighted!) by                                                */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*****************************************************************************/

/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

//#ifdef SINGLE
//#define REAL float
//#else /* not SINGLE */
//#define REAL double
//#endif /* not SINGLE */

#include <stdio.h>
#include <stdlib.h>
//#include "triangle.h"
#include "triangle.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <chrono>


// https://programmersought.com/article/46492593148/

//void report(struct triangulateio* io,int markers,int reporttriangles,int reportneighbors,int reportsegments,
//            int reportedges,int reportnorms)
//{
//    int i, j;
//
//    for (i = 0; i < io->numberofpoints; i++) {
//        printf("Point %4d:", i);
//        for (j = 0; j < 2; j++) {
//            printf("  %.6g", io->pointlist[i * 2 + j]);
//        }
//        if (io->numberofpointattributes > 0) {
//            printf("   attributes");
//        }
//        for (j = 0; j < io->numberofpointattributes; j++) {
//            printf("  %.6g",
//                   io->pointattributelist[i * io->numberofpointattributes + j]);
//        }
//        if (markers) {
//            printf("   marker %d\n", io->pointmarkerlist[i]);
//        } else {
//            printf("\n");
//        }
//    }
//    printf("\n");
//
//    if (reporttriangles || reportneighbors) {
//        for (i = 0; i < io->numberoftriangles; i++) {
//            if (reporttriangles) {
//                printf("Triangle %4d points:", i);
//                for (j = 0; j < io->numberofcorners; j++) {
//                    printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
//                }
////                if (io->numberoftriangleattributes > 0) {
////                    printf("   attributes");
////                }
////                for (j = 0; j < io->numberoftriangleattributes; j++) {
////                    printf("  %.6g", io->triangleattributelist[i *
////                                                               io->numberoftriangleattributes + j]);
////                }
//                printf("\n");
//            }
////            if (reportneighbors) {
////                printf("Triangle %4d neighbors:", i);
////                for (j = 0; j < 3; j++) {
////                    printf("  %4d", io->neighborlist[i * 3 + j]);
////                }
////                printf("\n");
////            }
//        }
//        printf("\n");
//    }
//
//    if (reportsegments) {
//        for (i = 0; i < io->numberofsegments; i++) {
//            printf("Segment %4d points:", i);
//            for (j = 0; j < 2; j++) {
//                printf("  %4d", io->segmentlist[i * 2 + j]);
//            }
//            if (markers) {
//                printf("   marker %d\n", io->segmentmarkerlist[i]);
//            } else {
//                printf("\n");
//            }
//        }
//        printf("\n");
//    }
//
////    if (reportedges) {
////        for (i = 0; i < io->numberofedges; i++) {
////            printf("Edge %4d points:", i);
////            for (j = 0; j < 2; j++) {
////                printf("  %4d", io->edgelist[i * 2 + j]);
////            }
////            if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
////                for (j = 0; j < 2; j++) {
////                    printf("  %.6g", io->normlist[i * 2 + j]);
////                }
////            }
////            if (markers) {
////                printf("   marker %d\n", io->edgemarkerlist[i]);
////            } else {
////                printf("\n");
////            }
////        }
////        printf("\n");
////    }
//}

std::vector<Eigen::Matrix<double,2,1>> readPolyVector(const std::string& fileName)
{
    std::vector<Eigen::Matrix<double,2,1>> temp;
    std::ifstream file ( fileName.c_str() , std::ifstream::in );
    if(file.is_open())
    {
        
        std::string line;
        double x,y;
        
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            
            ss >> x >>y;
            temp.push_back((Eigen::Matrix<double,2,1>()<<x,y).finished());
//            std::cout<<x<<","<<y<<std::endl;
            
            //poly[0] <<ClipperLib::DoublePoint(x,y);
        }
        
    }
    else
    {
        std::cout<<"CANNOT READ "+fileName<<std::endl;
    }
    
    
    
    return temp;
}

struct TriangleInterface //: public std::vector<Eigen::Matrix<double,2,1>>
{
    
//    TriangleInterface(const std::vector<Eigen::Matrix<double,2,1>>& pts) :
//    /* init */ std::vector<Eigen::Matrix<double,2,1>>(pts)
//    {
//
//    }
    
//    std::vector<Eigen::Matrix<double,2,1>>& points()
//    {
//        return *this;
//    }
//
//    const std::vector<Eigen::Matrix<double,2,1>>& points() const
//    {
//        return *this;
//    }
    
    void printVertices(const struct triangulateio* const io,const std::string& filename) const
    {
        std::ofstream verticesFile(filename);
        for (int i = 0; i < io->numberofpoints; i++)
        {
            verticesFile<<io->pointlist[i * 2 ]<<" "<<io->pointlist[i * 2 +1]<<"\n";

//            printf("Point %4d:", i);
//            for (j = 0; j < 2; j++)
//            {
//                vertices<<io->pointlist[i * 2 ]<<" "<<io->pointlist[i * 2 +1]<<"\n";
////                printf("  %.6g", io->pointlist[i * 2 + j]);
//            }
//            if (io->numberofpointattributes > 0) {
//                printf("   attributes");
//            }
//            for (j = 0; j < io->numberofpointattributes; j++) {
//                printf("  %.6g",
//                       io->pointattributelist[i * io->numberofpointattributes + j]);
//            }
//            if (markers) {
//                printf("   marker %d\n", io->pointmarkerlist[i]);
//            } else {
//                printf("\n");
//            }
        }
    }
    
    void printTriangles(const struct triangulateio* const io,const std::string& filename) const
    {
        std::ofstream trianglesFile(filename);
//        std::cout<<"#triangles="<<io->numberoftriangles<<std::endl;
        for (int i = 0; i < io->numberoftriangles; i++)
        {
////            if (reporttriangles) {
//                printf("Triangle %4d points:", i);
                for (int j = 0; j < io->numberofcorners; j++)
                {
//                    printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
                    trianglesFile<<io->trianglelist[i * io->numberofcorners + j]<<" ";
                }
            trianglesFile<<"\n";
////                if (io->numberoftriangleattributes > 0) {
////                    printf("   attributes");
////                }
////                for (j = 0; j < io->numberoftriangleattributes; j++) {
////                    printf("  %.6g", io->triangleattributelist[i *
////                                                               io->numberoftriangleattributes + j]);
////                }
//                printf("\n");
////            }
////            if (reportneighbors) {
////                printf("Triangle %4d neighbors:", i);
////                for (j = 0; j < 3; j++) {
////                    printf("  %4d", io->neighborlist[i * 3 + j]);
////                }
////                printf("\n");
////            }
        }

    }
    
    void mesh(const std::vector<Eigen::Matrix<double,2,1>>& boundaryPts,
              const std::vector<Eigen::Matrix<double,2,1>>& internalPts)
    {// https://userpages.umbc.edu/~rostamia/cbook/mesh.c
        // https://userpages.umbc.edu/~rostamia/cbook/triangle.h
        // https://people.sc.fsu.edu/~jburkardt/data/poly/poly.html
        std::cout<<"Creating triangular mesh"<<std::endl;
        struct triangulateio in,mid;// out, vorout;

//        numberofReservedpoints=100;
        
        in.numberofpoints = boundaryPts.size()+internalPts.size();
        in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
        for(int i=0;i<boundaryPts.size();++i)
        {
//            std::cout<<points()[p]<<std::endl;
            in.pointlist[2*i  ] = boundaryPts[i](0);
            in.pointlist[2*i+1] = boundaryPts[i](1);
        }
        for(int i=0;i<internalPts.size();++i)
        {
            const int i1(i+boundaryPts.size());
            //            std::cout<<points()[p]<<std::endl;
            in.pointlist[2*i1  ] = internalPts[i](0);
            in.pointlist[2*i1+1] = internalPts[i](1);
        }
//        in.pointlist[2*points().size()  ] = 0.5;
//        in.pointlist[2*points().size()+1] = 1.2;

//        in.pointlist[0] = 0.0;
//        in.pointlist[1] = 0.0;
//        in.pointlist[2] = 1.0;
//        in.pointlist[3] = 0.0;
//        in.pointlist[4] = 1.0;
//        in.pointlist[5] = 10.0;
//        in.pointlist[6] = 0.0;
//        in.pointlist[7] = 10.0;
        
        in.numberofpointattributes = 1;
        in.pointattributelist = (REAL *) malloc(in.numberofpoints *
                                                in.numberofpointattributes *
                                                sizeof(REAL));
        for(int i=0;i<in.numberofpoints;++i)
        {
            for(int j=0;j<in.numberofpointattributes;++j)
            {
                in.pointattributelist[i*in.numberofpointattributes+j] = 1.0;
            }
        }
//        in.pointattributelist[0] = 0.0;
//        in.pointattributelist[1] = 1.0;
//        in.pointattributelist[2] = 11.0;
//        in.pointattributelist[3] = 10.0;
        in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
        for(int i=0;i<in.numberofpoints;++i)
        {
            in.pointmarkerlist[i] = i;
        }
//        in.pointmarkerlist[0] = 0;
//        in.pointmarkerlist[1] = 2;
//        in.pointmarkerlist[2] = 0;
//        in.pointmarkerlist[3] = 0;
//        in.pointmarkerlist[4] = 2;
//        in.pointmarkerlist[5] = 2;
//        in.pointmarkerlist[6] = 0;

        
//                in.numberofsegments = 0;

        in.numberofsegments = boundaryPts.size();
        in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
//        in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
        in.segmentmarkerlist = (int *) NULL;

        for(int p=0;p<in.numberofsegments;++p)
        {
//                        std::cout<<p<<std::endl;
            in.segmentlist[2*p  ] = p;
            in.segmentlist[2*p+1] = (p<in.numberofsegments-1? p+1 : 0);
//            in.segmentmarkerlist[p]=p;
        }
        
        in.numberofholes = 0;
        
        
        float meshSize=10.0;
        in.numberofregions = 1;
        in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
//        in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
        in.regionlist[0] = 0.5;
        in.regionlist[1] = 0.5;
        in.regionlist[2] = 10.0;            /* Regional attribute (for whole mesh). */
        in.regionlist[3] = meshSize;          /* Area constraint that will not be used. */ // GP THIS CONTROLS TH MESH SIZE IN OUT???

    
//        printf("Input point set:\n\n");
//        report(&in, 1, 1, 1, 1, 1, 1);
        
        /* Make necessary initializations so that Triangle can return a */
        /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */
        
        mid.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
        /* Not needed if -N switch used or number of point attributes is zero: */
        mid.pointattributelist = (REAL *) NULL;
        mid.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
        mid.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
        /* Not needed if -E switch used or number of triangle attributes is zero: */
        mid.triangleattributelist = (REAL *) NULL;
        mid.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
        /* Needed only if segments are output (-p or -c) and -P not used: */
        mid.segmentlist = (int *) NULL;
        /* Needed only if segments are output (-p or -c) and -P and -B not used: */
        mid.segmentmarkerlist = (int *) NULL;
        mid.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
        mid.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */
//        mid.regionlist = (REAL *) NULL;
        mid.numberofregions = 0;
//        mid.regionlist = (REAL *) malloc(mid.numberofregions * 4 * sizeof(REAL));
//        mid.regionlist[0] = 0.5;
//        mid.regionlist[1] = 0.5;
//        mid.regionlist[2] = 10.0;            /* Regional attribute (for whole mesh). */
//        mid.regionlist[3] = 10.0;          /* Area constraint that will not be used. */
        
        
//        vorout.pointlist = (REAL *) NULL;        /* Needed only if -v switch used. */
//        /* Needed only if -v switch used and number of attributes is not zero: */
//        vorout.pointattributelist = (REAL *) NULL;
//        vorout.edgelist = (int *) NULL;          /* Needed only if -v switch used. */
//        vorout.normlist = (REAL *) NULL;         /* Needed only if -v switch used. */
        
        /* Triangulate the points.  Switches are chosen to read and write a  */
        /*   PSLG (p), preserve the convex hull (c), number everything from  */
        /*   zero (z), assign a regional attribute to each element (A), and  */
        /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
        /*   neighbor list (n).                                              */
        
//        triangulate("pzAevn", &in, &mid, &vorout); // G.P. use D for Delaunay triangulation, q to enforce angle quality
        std::cout<<"Meshing plane..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        triangulate("pazq", &in, &mid, (struct triangulateio *) NULL);
        std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;

        
                printVertices(&mid,"vertices_mid.txt");
                printTriangles(&mid,"triangles_mid.txt");

        
//        printf("Initial triangulation:\n\n");
//        report(&mid, 1, 1, 1, 1, 1, 1);
//        printf("Initial Voronoi diagram:\n\n");
//        report(&vorout, 0, 0, 0, 0, 1, 1);
        
        /* Attach area constraints to the triangles in preparation for */
        /*   refining the triangulation.                               */
        
        /* Needed only if -r and -a switches used: */
//        mid.trianglearealist = (REAL *) malloc(mid.numberoftriangles * sizeof(REAL));
//        float maxArea=10.0;
//        mid.trianglearealist[0] = maxArea;
//        mid.trianglearealist[1] = maxArea;
//        mid.trianglearealist[3] = maxArea;

        
        /* Make necessary initializations so that Triangle can return a */
        /*   triangulation in `out'.                                    */
     
//        out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
//        /* Not needed if -N switch used or number of attributes is zero: */
//        out.pointattributelist = (REAL *) NULL;
//        out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
//        /* Not needed if -E switch used or number of triangle attributes is zero: */
//        out.triangleattributelist = (REAL *) NULL;
//
//        /* Refine the triangulation according to the attached */
//        /*   triangle area constraints.                       */
//
//        triangulate("pazqBP", &mid, &out, (struct triangulateio *) NULL);
////        triangulate("pz", &mid, &out, (struct triangulateio *) NULL);
//
//        printf("Refined triangulation:\n\n");
////        report(&out, 0, 1, 0, 0, 0, 0);
//
//
//        printVertices(&out,"vertices_out.txt");
//        printTriangles(&out,"triangles_out.txt");

        
        /* Free all allocated arrays, including those allocated by Triangle. */
        
//        printf("Freeing in :\n\n");
        free(in.pointlist);
        free(in.pointattributelist);
        free(in.pointmarkerlist);
//        printf("Freeing mid :\n\n");
        free(in.regionlist);
        free(mid.pointlist);
        free(mid.pointattributelist);
        free(mid.pointmarkerlist);
        free(mid.trianglelist);
        free(mid.triangleattributelist);
//        free(mid.trianglearealist);
        free(mid.neighborlist);
        free(mid.segmentlist);
        free(mid.segmentmarkerlist);
        free(mid.edgelist);
        free(mid.edgemarkerlist);
//        printf("Freeing varout :\n\n");

//        free(vorout.pointlist);
//        free(vorout.pointattributelist);
//        free(vorout.edgelist);
//        free(vorout.normlist);
        
//        printf("Freeing out :\n\n");

//        free(out.pointlist);
//        free(out.pointattributelist);
//        free(out.trianglelist);
//        free(out.triangleattributelist);
    
    }
    
    
};

/*****************************************************************************/
/*                                                                           */
/*  report()   Print the input or output.                                    */
/*                                                                           */
/*****************************************************************************/

//void report(struct triangulateio* io,int markers,int reporttriangles,int reportneighbors,int reportsegments,
//            int reportedges,int reportnorms)
//{
//  int i, j;
//
//  for (i = 0; i < io->numberofpoints; i++) {
//    printf("Point %4d:", i);
//    for (j = 0; j < 2; j++) {
//      printf("  %.6g", io->pointlist[i * 2 + j]);
//    }
//    if (io->numberofpointattributes > 0) {
//      printf("   attributes");
//    }
//    for (j = 0; j < io->numberofpointattributes; j++) {
//      printf("  %.6g",
//             io->pointattributelist[i * io->numberofpointattributes + j]);
//    }
//    if (markers) {
//      printf("   marker %d\n", io->pointmarkerlist[i]);
//    } else {
//      printf("\n");
//    }
//  }
//  printf("\n");
//
//  if (reporttriangles || reportneighbors) {
//    for (i = 0; i < io->numberoftriangles; i++) {
//      if (reporttriangles) {
//        printf("Triangle %4d points:", i);
//        for (j = 0; j < io->numberofcorners; j++) {
//          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
//        }
//        if (io->numberoftriangleattributes > 0) {
//          printf("   attributes");
//        }
//        for (j = 0; j < io->numberoftriangleattributes; j++) {
//          printf("  %.6g", io->triangleattributelist[i *
//                                         io->numberoftriangleattributes + j]);
//        }
//        printf("\n");
//      }
//      if (reportneighbors) {
//        printf("Triangle %4d neighbors:", i);
//        for (j = 0; j < 3; j++) {
//          printf("  %4d", io->neighborlist[i * 3 + j]);
//        }
//        printf("\n");
//      }
//    }
//    printf("\n");
//  }
//
//  if (reportsegments) {
//    for (i = 0; i < io->numberofsegments; i++) {
//      printf("Segment %4d points:", i);
//      for (j = 0; j < 2; j++) {
//        printf("  %4d", io->segmentlist[i * 2 + j]);
//      }
//      if (markers) {
//        printf("   marker %d\n", io->segmentmarkerlist[i]);
//      } else {
//        printf("\n");
//      }
//    }
//    printf("\n");
//  }
//
//  if (reportedges) {
//    for (i = 0; i < io->numberofedges; i++) {
//      printf("Edge %4d points:", i);
//      for (j = 0; j < 2; j++) {
//        printf("  %4d", io->edgelist[i * 2 + j]);
//      }
//      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
//        for (j = 0; j < 2; j++) {
//          printf("  %.6g", io->normlist[i * 2 + j]);
//        }
//      }
//      if (markers) {
//        printf("   marker %d\n", io->edgemarkerlist[i]);
//      } else {
//        printf("\n");
//      }
//    }
//    printf("\n");
//  }
//}

/*****************************************************************************/
/*                                                                           */
/*  main()   Create and refine a mesh.                                       */
/*                                                                           */
/*****************************************************************************/

int main()
{
    
    std::vector<Eigen::Matrix<double,2,1>> boundaryPoints=readPolyVector("boundaryPoints.txt");
    std::vector<Eigen::Matrix<double,2,1>> internalPoints=readPolyVector("internalPoints.txt");

    TriangleInterface ti;
    ti.mesh(boundaryPoints,internalPoints);
    
    
//  struct triangulateio in, mid, out, vorout;
//
//  /* Define input points. */
//
//  in.numberofpoints = 4;
//  in.numberofpointattributes = 1;
//  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
//  in.pointlist[0] = 0.0;
//  in.pointlist[1] = 0.0;
//  in.pointlist[2] = 1.0;
//  in.pointlist[3] = 0.0;
//  in.pointlist[4] = 1.0;
//  in.pointlist[5] = 10.0;
//  in.pointlist[6] = 0.0;
//  in.pointlist[7] = 10.0;
//  in.pointattributelist = (REAL *) malloc(in.numberofpoints *
//                                          in.numberofpointattributes *
//                                          sizeof(REAL));
//  in.pointattributelist[0] = 0.0;
//  in.pointattributelist[1] = 1.0;
//  in.pointattributelist[2] = 11.0;
//  in.pointattributelist[3] = 10.0;
//  in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
//  in.pointmarkerlist[0] = 0;
//  in.pointmarkerlist[1] = 2;
//  in.pointmarkerlist[2] = 0;
//  in.pointmarkerlist[3] = 0;
//
//  in.numberofsegments = 0;
//  in.numberofholes = 0;
//  in.numberofregions = 1;
//  in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
//  in.regionlist[0] = 0.5;
//  in.regionlist[1] = 5.0;
//  in.regionlist[2] = 7.0;            /* Regional attribute (for whole mesh). */
//  in.regionlist[3] = 0.1;          /* Area constraint that will not be used. */
//
//  printf("Input point set:\n\n");
//  report(&in, 1, 0, 0, 0, 0, 0);
//
//  /* Make necessary initializations so that Triangle can return a */
//  /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */
//
//  mid.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
//  /* Not needed if -N switch used or number of point attributes is zero: */
//  mid.pointattributelist = (REAL *) NULL;
//  mid.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
//  mid.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
//  /* Not needed if -E switch used or number of triangle attributes is zero: */
//  mid.triangleattributelist = (REAL *) NULL;
//  mid.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
//  /* Needed only if segments are output (-p or -c) and -P not used: */
//  mid.segmentlist = (int *) NULL;
//  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
//  mid.segmentmarkerlist = (int *) NULL;
//  mid.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
//  mid.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */
//
//  vorout.pointlist = (REAL *) NULL;        /* Needed only if -v switch used. */
//  /* Needed only if -v switch used and number of attributes is not zero: */
//  vorout.pointattributelist = (REAL *) NULL;
//  vorout.edgelist = (int *) NULL;          /* Needed only if -v switch used. */
//  vorout.normlist = (REAL *) NULL;         /* Needed only if -v switch used. */
//
//  /* Triangulate the points.  Switches are chosen to read and write a  */
//  /*   PSLG (p), preserve the convex hull (c), number everything from  */
//  /*   zero (z), assign a regional attribute to each element (A), and  */
//  /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
//  /*   neighbor list (n).                                              */
//
//  triangulate("pczAevn", &in, &mid, &vorout);
//
//  printf("Initial triangulation:\n\n");
//  report(&mid, 1, 1, 1, 1, 1, 0);
//  printf("Initial Voronoi diagram:\n\n");
//  report(&vorout, 0, 0, 0, 0, 1, 1);
//
//  /* Attach area constraints to the triangles in preparation for */
//  /*   refining the triangulation.                               */
//
//  /* Needed only if -r and -a switches used: */
//  mid.trianglearealist = (REAL *) malloc(mid.numberoftriangles * sizeof(REAL));
//  mid.trianglearealist[0] = 3.0;
//  mid.trianglearealist[1] = 1.0;
//
//  /* Make necessary initializations so that Triangle can return a */
//  /*   triangulation in `out'.                                    */
//
//  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
//  /* Not needed if -N switch used or number of attributes is zero: */
//  out.pointattributelist = (REAL *) NULL;
//  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
//  /* Not needed if -E switch used or number of triangle attributes is zero: */
//  out.triangleattributelist = (REAL *) NULL;
//
//  /* Refine the triangulation according to the attached */
//  /*   triangle area constraints.                       */
//
//  triangulate("prazBP", &mid, &out, (struct triangulateio *) NULL);
//
//  printf("Refined triangulation:\n\n");
//  report(&out, 0, 1, 0, 0, 0, 0);
//
//  /* Free all allocated arrays, including those allocated by Triangle. */
//
//  free(in.pointlist);
//  free(in.pointattributelist);
//  free(in.pointmarkerlist);
//  free(in.regionlist);
//  free(mid.pointlist);
//  free(mid.pointattributelist);
//  free(mid.pointmarkerlist);
//  free(mid.trianglelist);
//  free(mid.triangleattributelist);
//  free(mid.trianglearealist);
//  free(mid.neighborlist);
//  free(mid.segmentlist);
//  free(mid.segmentmarkerlist);
//  free(mid.edgelist);
//  free(mid.edgemarkerlist);
//  free(vorout.pointlist);
//  free(vorout.pointattributelist);
//  free(vorout.edgelist);
//  free(vorout.normlist);
//  free(out.pointlist);
//  free(out.pointattributelist);
//  free(out.trianglelist);
//  free(out.triangleattributelist);

  return 0;
}

//
// Created by milinda on 9/29/16.
//
/**
 *
 * @author Milinda Fernando
 * School of Computing, University of Utah.
 * @breif Contains the basic example to generate the octrees for the given black holes.
 *
 * */

#include "mpi.h"
#include "TreeNode.h"
#include "sfcSort.h"
#include "octUtils.h"
#include "interpolate.h"
#include "visualize.h"


int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int rank, npes;
    MPI_Comm GLOBAL_COMM = MPI_COMM_WORLD;
    MPI_Comm_rank(GLOBAL_COMM, &rank);
    MPI_Comm_size(GLOBAL_COMM, &npes);

    //std::cout<<"rank: "<<rank<<" npes: "<<npes<<std::endl;

    unsigned int options = 0;
    auto t1=std::chrono::high_resolution_clock::now();
    double otGeneration=0;
    double otBalancing=0;


    if (argc < 8) {
        if (!rank)
            std::cout << "Usage :" << argv[0] << "wtol(max) iterations (wtol(min)=wtol(max)/10^iterations)" << " dim " << " maxDepth " << " tsortTol" << " options( 1-Remove duplicates 2- construct octree 6-balanced octree) SF_K(=128) VTKFileName=(balOct)"<< std::endl;
    }

    double wavelet_tol=1e-3;

    int iterations=1;
    unsigned int sf_k=128;
    char vtkName[256];
    wavelet_tol=atof(argv[1]);
    iterations=atoi(argv[2]);


    unsigned int dim=atoi(argv[3]);
    unsigned int maxDepth=atoi(argv[4]);
    double tol=0.001;
    unsigned int distribution=0;
    if(argc>5)
        tol=atof(argv[5]);
    if(argc>6)
        options=atoi(argv[6]);
    if(argc>7)
        sf_k=atoi(argv[7]);



    _InitializeHcurve(dim);
    if(!rank) std::cout << "Initialized Hcurves for dimention "<<dim << " wavelet tol (max): "<<wavelet_tol<<" iterations: "<<iterations<<std::endl;



    DendroIntL localSz, totalSz, totalBalSz;
    std::vector<ot::TreeNode> tmpNodes;

    unsigned int inital_lev=std::ceil(binOp::fastLog2(npes)/(double)3);

    createRegularOctree(tmpNodes,inital_lev,dim,maxDepth,MPI_COMM_WORLD);
    //treeNodesTovtk(tmpNodes,rank,"initalOctree");
    std::vector<ot::TreeNode> refined;
    std::vector<ot::TreeNode> children;
    std::vector<ot::TreeNode> tmpSorted;
    std::vector<ot::TreeNode> balOct;

    for(unsigned int iter=0;iter<iterations;iter++) {


        if(!rank) std::cout<<" iter: "<<iter<<" wtol: "<<wavelet_tol<<std::endl;
        sprintf(vtkName,"%s_tol_%f",argv[8],wavelet_tol);
        refined.clear();
        children.clear();
        tmpSorted.clear();
        balOct.clear();

        ot::TreeNode tmp;
        t1=std::chrono::high_resolution_clock::now();
        for (unsigned int i = 0; i < tmpNodes.size(); i++) {
            if (isSplit_linear(tmpNodes[i], wavelet_tol)) {
                tmpNodes[i].addChildren(children);
                refined.insert(refined.end(), children.begin(), children.end());
                tmpNodes.insert(tmpNodes.end(), children.begin(), children.end());
                children.clear();
            } else {
                refined.push_back(tmpNodes[i]);
            }
            /*if(!rank)
                std::cout<<"tmpNode Size: "<<tmpNodes.size()<<std::endl;*/
        }

        if(!rank) std::cout<<"refinement complete"<<std::endl;
        tmpNodes.clear();
        ot::TreeNode root = ot::TreeNode(0, 0, 0, 0, dim, maxDepth);

        SFC::parSort::SFC_treeSort(refined, tmpSorted, tmpSorted, tmpSorted, tol, maxDepth, root, 0, 1, options, sf_k,
                                   MPI_COMM_WORLD);
        otGeneration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - t1).count();

        //std::swap(tmpNodes,tmp);
        refined.clear();

        if(!rank) std::cout<<"octree construction completed"<<std::endl;
        localSz = tmpSorted.size();
        par::Mpi_Reduce(&localSz, &totalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);

        /* if(!rank)
         {
              std::cout<<"=============== Properties ============================="<<std::endl;
              std::cout<<" Wavelet Tolerance: "<<wavelet_tol<<std::endl;
              std::cout<<" Total octree(before 2:1 balance) elements: "<<totalSz<<std::endl;
              std::cout<<"========================================================"<<std::endl;
         }*/



        t1 = std::chrono::high_resolution_clock::now();
        SFC::parSort::SFC_treeSort(tmpSorted, balOct, balOct, balOct, tol, maxDepth, root, 0, 1, 4, sf_k, MPI_COMM_WORLD);
        otBalancing = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - t1).count();


        localSz = balOct.size();
        par::Mpi_Reduce(&localSz, &totalBalSz, 1, MPI_SUM, 0, MPI_COMM_WORLD);
        //std::swap(tmpSorted, balOct);
        //balOct.clear();
        std::swap(tmpNodes,tmpSorted);
        tmpSorted.clear();
        if(!rank) std::cout<<"octree 2:1 balancing completed"<<std::endl;
        double otGen_g[3]; // min mean max
        double otBal_g[3]; // min mean max

        par::Mpi_Reduce(&otGeneration, otGen_g, 1, MPI_MIN, 0, MPI_COMM_WORLD);
        par::Mpi_Reduce(&otGeneration, otGen_g + 1, 1, MPI_SUM, 0, MPI_COMM_WORLD);
        par::Mpi_Reduce(&otGeneration, otGen_g + 2, 1, MPI_MAX, 0, MPI_COMM_WORLD);
        otGen_g[1] = otGen_g[1] / (double) npes;

        par::Mpi_Reduce(&otBalancing, otBal_g, 1, MPI_MIN, 0, MPI_COMM_WORLD);
        par::Mpi_Reduce(&otBalancing, otBal_g + 1, 1, MPI_SUM, 0, MPI_COMM_WORLD);
        par::Mpi_Reduce(&otBalancing, otBal_g + 2, 1, MPI_MAX, 0, MPI_COMM_WORLD);
        otBal_g[1] = otBal_g[1] / (double) npes;

        otGen_g[0] = otGen_g[0] / (double) 1000;
        otGen_g[1] = otGen_g[1] / (double) 1000;
        otGen_g[2] = otGen_g[2] / (double) 1000;

        otBal_g[0] = otBal_g[0] / (double) 1000;
        otBal_g[1] = otBal_g[1] / (double) 1000;
        otBal_g[2] = otBal_g[2] / (double) 1000;


        if (!rank) {
            std::cout
                    << "npes\twavelet_tol\tmaxDepth\toctSz\tbalOctSz\toctGen_min\toctGen_mean\toctGen_max\toctBal_min\toctBal_mean\toctBal_max"
                    << std::endl;
            std::cout << npes << "\t" << wavelet_tol << "\t" << maxDepth << "\t" << totalSz << "\t" << totalBalSz
                      << "\t" << otGen_g[0] << "\t" << otGen_g[1] << "\t" << otGen_g[2] << "\t" << otBal_g[0] << "\t"
                      << otBal_g[1] << "\t" << otBal_g[2] << std::endl;
        }


        double *u = new double[NUM_FUNC_VALUES * balOct.size() * NUM_CHILDREN];
        Point p[NUM_CHILDREN];
        unsigned int myX, myY, myZ, mySz;

        for (unsigned int i = 0; i < balOct.size(); i++) {
            myX = balOct[i].getX();
            myY = balOct[i].getY();
            myZ = balOct[i].getZ();

            mySz = 1u << (m_uiMaxDepth - balOct[i].getLevel());

            if (m_uiDim == 2) {
                p[0] = Point(myX, myY, myZ);
                p[1] = Point(myX + mySz, myY, myZ);
                p[2] = Point(myX + mySz, myY + mySz, myZ);
                p[3] = Point(myX, myY + mySz, myZ);
            } else if (m_uiDim == 3) {
                p[0] = Point(myX, myY, myZ);
                p[1] = Point(myX + mySz, myY, myZ);
                p[2] = Point(myX + mySz, myY + mySz, myZ);
                p[3] = Point(myX, myY + mySz, myZ);

                p[4] = Point(myX, myY, myZ + mySz);
                p[5] = Point(myX + mySz, myY, myZ + mySz);
                p[6] = Point(myX + mySz, myY + mySz, myZ + mySz);
                p[7] = Point(myX, myY + mySz, myZ + mySz);

            }
            for (unsigned int j = 0; j < NUM_CHILDREN; j++)
                init_data_puncture(p[j], (u + i * (NUM_CHILDREN * NUM_FUNC_VALUES) + j * NUM_FUNC_VALUES));
        }

        //treeNodesTovtk(tmpSorted,rank,"refined");
        //octree2VTK(tmpSorted,rank,u,vtkName);
        //if(!rank) std::cout<<"balOct was written to : "<<vtkName<<std::endl;
        //octree2VTK(tmpSorted,rank,u,"refinedOt");
        balOct.clear();
        delete[] u;
        wavelet_tol=wavelet_tol/10.0;
    }


    MPI_Finalize();

}
//
// Created by milinda on 10/6/16.
//

#include "visualize.h"

/**
 * @author: Milinda Fernando
 * School of Computing , University of Utah.
 *
 * @brief: Contains the functions related to visulizing octress, and mesh data in vtk format.
 *
 * */

/**
 * @author Milinda Fernando
 * @breif Writes the octree to the vtk file format with the function values. (currenlty all the 24 functions. )
 * @param[in] octree -> octree which needs to be written
 * @param[in] rank -> mpi rank
 * @param[in] function values corresponding to each octant.
 */

void octree2VTK(const std::vector<ot::TreeNode>& nodes, unsigned int rank, const double * u,std::string vtk_file_name)
{

    if (!rank) std::cout << "writing mesh to VTK file: " << vtk_file_name << std::endl;
    std::ostringstream convert;

#ifdef  HILBERT_ORDERING
    convert << vtk_file_name << "_H_" <<rank << ".vtk";
#else
    convert << vtk_file_name << "_M_" << mpi_rank << ".vtk";
#endif

    //convert << vtk_file_name << "_" << mpi_rank << ".vtk";
    vtk_file_name = convert.str();

    std::ofstream myfile;
    myfile.open(vtk_file_name.c_str());

    myfile << "# vtk DataFile Version 2.0" << std::endl;
    myfile << "DENDRO OCTREES" << std::endl;
    myfile << "ASCII" << std::endl;
    myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    int dim = 3;//nodes[0].getDim();

    int unit_points = 1 << dim;
    int num_verticies = nodes.size() * (unit_points);
    int num_cells = nodes.size();

    //std::cout<<"FileName :"<<vtk_file_name<<"\t DIM:"<<dim<<"\t number of points:"<<num_verticies<<"\t number of cells:"<<num_cells<<std::endl;


    myfile << "POINTS " << num_verticies << " float" << std::endl;

    if (dim == 2) {

    } else if (dim == 3) {
        unsigned int len;
        unsigned int xl, yl, zl;
        int num_data_field = 2;
        /* if (hsorted) {
             num_data_field++;
             std::sort(nodes.begin(), nodes.end());
         }*/

        int num_cells_elements = num_cells * unit_points + num_cells;

        for (int i = 0; i < nodes.size(); i++) {
            //std::cout<<nodes[i]<<std::endl;
            len = 1 << (nodes[i].getMaxDepth() - nodes[i].getLevel());
            xl = nodes[i].getX();
            yl = nodes[i].getY();
            zl = nodes[i].getZ();

            myfile << xl << " " << yl << " " << zl << std::endl;
            myfile << (xl + len) << " " << yl << " " << zl << std::endl;
            myfile << (xl + len) << " " << (yl + len) << " " << zl << std::endl;
            myfile << xl << " " << (yl + len) << " " << zl << std::endl;

            myfile << xl << " " << yl << " " << (zl + len) << std::endl;
            myfile << (xl + len) << " " << yl << " " << (zl + len) << std::endl;
            myfile << (xl + len) << " " << (yl + len) << " " << (zl + len) << std::endl;
            myfile << xl << " " << (yl + len) << " " << (zl + len) << std::endl;

        }

        myfile << "CELLS " << nodes.size() << " " << num_cells_elements << std::endl;

        for (int i = 0; i < num_cells; i++) {
            myfile << unit_points << " ";
            for (int j = 0; j < unit_points; j++) {
                myfile << (i * unit_points + j) << " ";
            }
            myfile << std::endl;
        }

        myfile << "CELL_TYPES " << num_cells << std::endl;
        for (int i = 0; i < num_cells; i++) {
            myfile << VTK_HEXAHEDRON << std::endl;
        }

        myfile<<std::endl;

        myfile<< "POINT_DATA "<<num_verticies<<std::endl;


        myfile <<"SCALARS U_ALPHA FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_ALPHA] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_CHI FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_CHI] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_TRK FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_TRK] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"VECTORS U_SHIFT FLOAT"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_SHIFTX]<<" "<<u[(j) * (NUM_FUNC_VALUES) + U_SHIFTY]<<" "<<u[(j) * (NUM_FUNC_VALUES) + U_SHIFTZ]<< std::endl;
              /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"VECTORS U_GMAT FLOAT"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_GAMTX]<<" "<<u[(j) * (NUM_FUNC_VALUES) + U_GAMTY]<<" "<<u[(j) * (NUM_FUNC_VALUES) + U_GAMTZ]<< std::endl;
            /*if(!rank)
              std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;


        myfile <<"VECTORS U_GB FLOAT"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_GBX]<<" "<<u[(j) * (NUM_FUNC_VALUES) + U_GBY]<<" "<<u[(j) * (NUM_FUNC_VALUES) + U_GBZ]<< std::endl;
            /*if(!rank)
              std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;


        myfile <<"SCALARS U_GTXX FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_GTXX] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_GTXY FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_GTXY] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_GTXZ FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_GTXZ] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;


        myfile <<"SCALARS U_GTYY FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_GTYY] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_GTYZ FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_GTYZ] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_GTZZ FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_GTZZ] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;


        myfile <<"SCALARS U_ATXX FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_ATXX] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_ATXY FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_ATXY] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_ATXZ FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_ATXZ] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;


        myfile <<"SCALARS U_ATYY FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_ATYY] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_ATYZ FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_ATYZ] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile <<"SCALARS U_ATZZ FLOAT"<<std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for(unsigned int j=0;j<num_verticies;j++) {
            myfile << u[(j) * (NUM_FUNC_VALUES) + U_ATZZ] << std::endl;
            /*if(!rank)
                std::cout<<"index: "<<((j) * (NUM_FUNC_VALUES) + i)<<" value: "<<u[(j) * (NUM_FUNC_VALUES) + i]<<std::endl;*/
        }
        myfile<<std::endl;

        myfile<<"CELL_DATA "<<num_cells<<std::endl;
        //myfile<<"POINT_DATA "<<(num_cells*unit_points)<<std::endl;



        //myfile << "FIELD OCTREE_DATA " << num_data_field << std::endl;

        //myfile << "cell_level 1 " << num_cells << " int" << std::endl;
        myfile << "SCALARS cell_level FLOAT"<< std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;

        for (int i = 0; i < nodes.size(); i++)
            myfile << nodes[i].getLevel() <<std::endl;

        myfile << std::endl;
        //myfile << std::endl;

        //myfile << "mpi_rank 1 " << num_cells << " int" << std::endl;
        myfile << "SCALARS mpi_rank FLOAT"<< std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for (int i = 0; i < nodes.size(); i++)
            myfile <<rank << std::endl;

        myfile << std::endl;






    }

    myfile.close();



}

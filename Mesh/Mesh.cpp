#include "Mesh.h"



Mesh::Mesh(std::string filename)
{
    std::cout<<"Starting mesh initialization. \n";

    std::ifstream file(filename);

    if(!file)
    {
        std::cout<<"ERROR: WRONG FILE NAME! \n";
    }
    
    unsigned int Nx, Ny, N_cells, N_faces;

    file >> Nx >> Ny >> N_cells >> N_faces;
    std::cout<< N_faces <<"\n";

    Cells.reserve(N_cells);
    Faces.reserve(N_faces);

    for(unsigned int k = 0; k < N_cells; ++k )
    {
        unsigned int id, i, j;
        unsigned int faces[4];

        file >> id >> i >> j >> faces[0] >> faces[1] >> faces[2] >> faces[3];

        Cells.emplace_back(id, i, j, faces);
    }

    for(unsigned int k = 0; k < N_faces; ++k )
    {
        unsigned int id, left_id, right_id;
        bool isVertiacal;
        int type;
        file >> id >> left_id >> right_id >> isVertiacal >> type;
        if(isVertiacal)
            VertToHoriz = k;

        Faces.emplace_back(id, left_id, right_id, isVertiacal, type);
    }


    std::cout<<"Mesh initialization finished. \n";

    file.close();


}

//Необходима реализация
void Mesh::SetInitialCondidions()
{
    unsigned int discontinuity = 600; //100
    double a1L, ro1L, u1L, v1L, P1L, ro2L, u2L, v2L, P2L,
           a1R, ro1R, u1R, v1R, P1R, ro2R, u2R, v2R, P2R;
    double WL[7];
//Сделать чтение из файла?
    a1L = 1e-6;
    ro1L = 2;
    u1L = 0;
    v1L = 0;
    P1L = 100000000.0;
    ro2L = 1150;
    u2L = 0;
    v2L = 0;
    P2L = 100000000.0;

    a1R = 1-1e-6;
    ro1R = 2;
    u1R = 0;
    v1R = 0;
    P1R = 100000.0;
    ro2R = 1150;
    u2R = 0;
    v2R = 0;
    P2R = 100000.0;


    StateW stateL(a1L, ro1L, u1L, v1L, P1L, ro2L, u2L, v2L, P2L);
    StateW stateR(a1R, ro1R, u1R, v1R, P1R, ro2R, u2R, v2R, P2R);
    for(
        auto& cell : Cells)
    {
        if(cell.i < discontinuity )
            cell.W = stateL;
        else 
            cell.W = stateR;
    }

}
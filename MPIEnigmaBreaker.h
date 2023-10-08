//
// Created by Mateusz Chyzy on 28/10/2022.
//

#ifndef UNTITLED_MPIENIGMABREAKER_H
#define UNTITLED_MPIENIGMABREAKER_H


#include <mpi.h>
#include "Consts.h"
#include <vector>
#include "EnigmaBreaker.h"

class MPIEnigmaBreaker : public EnigmaBreaker {
    const int ROOT = 0;
    const int SIZE_BUFFER = 500000;
    const int TAG_SOLUTION_FOUND = 100;
    const int TAG_ROOT = 0;
    int flag = 0;
    int rank;
    int size;
    uint *expectedMessage;
    uint expectedLength = 0;
    uint *messageToDecode;
    uint messageLength = 0;
public:
    MPIEnigmaBreaker(Enigma *enigma, MessageComparator *comparator);

    void crackMessage();
    virtual void setSampleToFind( uint *expected, uint expectedLenght );
    virtual void setMessageToDecode( uint *message, uint messageLength );
    void calculateRange(vector<vector<uint>> &tmp);
    void findRotorsPositions();

    int MPI_Comm_rank(MPI_Comm comm, int *rank);
    int MPI_Barrier(MPI_Comm comm);

    bool solutionFound( uint *rotorSettingsProposal );

    virtual void getResult( uint *rotorPositions ) override;
    virtual ~MPIEnigmaBreaker();
};


#endif //UNTITLED_MPIENIGMABREAKER_H
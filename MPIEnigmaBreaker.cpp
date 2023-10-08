//
// Created by Mateusz Chyzy on 01/10/2023.
//
#include <vector>
#include "MPIEnigmaBreaker.h"

using namespace std;

MPIEnigmaBreaker::MPIEnigmaBreaker(Enigma *enigma, MessageComparator *comparator) : EnigmaBreaker(enigma, comparator) {
   MPIEnigmaBreaker::MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
   MPI_Comm_size(MPI_COMM_WORLD, &this->size);
}

void MPIEnigmaBreaker::crackMessage() {
   char buffer[SIZE_BUFFER];
   int position;

   if (rank == ROOT) {
       position = 0;

       MPI_Pack(&messageLength, 1, MPI_UNSIGNED, buffer, SIZE_BUFFER, &position, MPI_COMM_WORLD);
       MPI_Pack(messageToDecode, this->messageLength, MPI_UNSIGNED, buffer, SIZE_BUFFER, &position, MPI_COMM_WORLD);
       MPI_Pack(&expectedLength, 1, MPI_UNSIGNED, buffer, SIZE_BUFFER, &position, MPI_COMM_WORLD);
       MPI_Pack(expectedMessage, this->expectedLength, MPI_UNSIGNED, buffer, SIZE_BUFFER, &position, MPI_COMM_WORLD);

       MPI_Bcast(buffer, SIZE_BUFFER, MPI_PACKED, ROOT, MPI_COMM_WORLD);

       findRotorsPositions();

   } else {
       MPI_Bcast(buffer, SIZE_BUFFER, MPI_PACKED, ROOT, MPI_COMM_WORLD);

       position = 0;
       MPI_Unpack(buffer, SIZE_BUFFER, &position, &this->messageLength, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
       this->messageToDecode = new uint[messageLength];

       MPI_Unpack(buffer, SIZE_BUFFER, &position, this->messageToDecode, this->messageLength, MPI_UNSIGNED,
                  MPI_COMM_WORLD);

       this->setMessageToDecode(this->messageToDecode, this->messageLength);

       MPI_Unpack(buffer, SIZE_BUFFER, &position, &this->expectedLength, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
       this->expectedMessage = new uint[expectedLength];

       MPI_Unpack(buffer, SIZE_BUFFER, &position, this->expectedMessage, this->expectedLength, MPI_UNSIGNED,
                  MPI_COMM_WORLD);
       EnigmaBreaker::setSampleToFind(this->expectedMessage, this->expectedLength);

       findRotorsPositions();

   }
   delete[] messageToDecode;
   delete[] expectedMessage;
}

bool MPIEnigmaBreaker::solutionFound(uint *rotorSettingsProposal) {
   for (uint rotor = 0; rotor < rotors; rotor++) {
       rotorPositions[rotor] = rotorSettingsProposal[rotor];
   }

   enigma->setRotorPositions(rotorPositions);
   uint *decodedMessage = new uint[messageLength];

   for (uint messagePosition = 0; messagePosition < messageLength; messagePosition++) {
       decodedMessage[messagePosition] = enigma->code(messageToDecode[messagePosition]);
   }
   bool result = comparator->messageDecoded(decodedMessage);

   delete[] decodedMessage;

   if (result) {
       MPI_Send(rotorPositions, rotors, MPI_UNSIGNED, ROOT, TAG_ROOT, MPI_COMM_WORLD);
   }

   return result;
}

void MPIEnigmaBreaker::setSampleToFind(uint *expected, uint expectedLenght) {
   this->expectedMessage = expected;
   this->expectedLength = expectedLenght;
   EnigmaBreaker::setSampleToFind(expected, expectedLenght);
}

void MPIEnigmaBreaker::setMessageToDecode(uint *message, uint messageLength) {
   this->messageToDecode = message;
   this->messageLength = messageLength;
   EnigmaBreaker::setMessageToDecode(message, messageLength);
}

MPIEnigmaBreaker::~MPIEnigmaBreaker() {
   delete[] rotorPositions;
   delete[] messageToDecode;
   delete[] expectedMessage;
}

void MPIEnigmaBreaker::getResult(uint *rotorPositions) {
   if (rank == ROOT) {
       MPI_Recv(this->rotorPositions, rotors, MPI_UNSIGNED, MPI_ANY_SOURCE, ROOT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       for (uint rotor = 0; rotor < rotors; rotor++) {
           rotorPositions[rotor] = this->rotorPositions[rotor];
       }
   }
}

int MPIEnigmaBreaker::MPI_Comm_rank(MPI_Comm comm, int *rank) {
   return ::MPI_Comm_rank(comm, rank);
}

void MPIEnigmaBreaker::findRotorsPositions() {
   uint rotorLargestSetting = enigma->getLargestRotorSetting();

   uint *rMax = new uint[MAX_ROTORS];
   for (uint rotor = 0; rotor < MAX_ROTORS; rotor++) {
       if (rotor < rotors)
           rMax[rotor] = rotorLargestSetting;
       else
           rMax[rotor] = 0;
   }

   MPI_Request req;
   MPI_Irecv(&this->flag, 1, MPI_INT, MPI_ANY_SOURCE, TAG_SOLUTION_FOUND, MPI_COMM_WORLD, &req);
   uint *r = new uint[ MAX_ROTORS ];

   for (r[0] = this->rank; r[0] <= rMax[0]; r[0]+=this->size) {
       for (r[1] = 0; r[1] <= rMax[1]; r[1]++) {
           for (r[2] = 0; r[2] <= rMax[2]; r[2]++) {
               for (r[3] = 0; r[3] <= rMax[3]; r[3]++) {
                   for (r[4] = 0; r[4] <= rMax[4]; r[4]++) {
                       for (r[5] = 0; r[5] <= rMax[5]; r[5]++) {
                           for (r[6] = 0; r[6] <= rMax[6]; r[6]++) {
                               for (r[7] = 0; r[7] <= rMax[7]; r[7]++) {
                                   for (r[8] = 0; r[8] <= rMax[8]; r[8]++) {
                                       for (r[9] = 0; r[9] <= rMax[9]; r[9]++) {
                                           if (!this->flag) {
                                               if (solutionFound(r)){
                                                   this->flag = 1;
                                                   for (size_t f = 0; f < size; f++) {
                                                       if (f == rank) {
                                                           continue;
                                                       } else {
                                                           MPI_Isend(&this->flag, 1, MPI_INT, f, TAG_SOLUTION_FOUND, MPI_COMM_WORLD, &req);
                                                       }
                                                   }
                                                   break;
                                               }
                                           }
                                           MPI_Test(&req, &this->flag, MPI_STATUS_IGNORE);
                                           if (this->flag) {
                                               goto EXIT_ALL_LOOPS;
                                           }
                                       }
                                   }
                               }
                           }
                       }
                   }
               }
           }
       }
   }
   EXIT_ALL_LOOPS:

   delete[] rMax;
   delete[] r;
}


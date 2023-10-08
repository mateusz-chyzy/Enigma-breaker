/*
 * Main.cpp
 */

#include <stdlib.h>
#include <iostream>
#include "mpi.h" 

#include"Consts.h"
#include"Enigma.h"
#include"Machinery.h"
#include"SimpleMachinery.h"
#include"MessageComparator.h"
#include"SimpleMessageComparator.h"
#include"EnigmaBreaker.h"
#include"SerialEnigmaBreaker.h"

const uint ROTORS = 6;
const uint LARGEST_ROTOR_SETTING = 10;
const uint MESSAGE_LENGTH = 1536;
const uint EXPECTED_MESSAGE_LENGTH = 64;
const uint MAX_VALUE_IN_MESSAGE = 256;

using namespace std;

void generateRandomMessage( uint *message, uint length ) {
	for ( uint i = 0; i < length; i++ )
		message[ i ] = random() % MAX_VALUE_IN_MESSAGE;
}

void copy( uint *source, uint sourceLength, uint* destination, uint destinationLength ) {
	uint shift = sourceLength - destinationLength;

	for ( uint pos = 0; pos < destinationLength; pos++ )
		destination[ pos ] = source[ pos + shift ];
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	Machinery *enigmaMachinery = new SimpleMachinery(ROTORS,LARGEST_ROTOR_SETTING);
	Enigma *enigma = new Enigma( enigmaMachinery );

	MessageComparator *comparator = new SimpleMessageComparator();
	EnigmaBreaker *breaker = new SerialEnigmaBreaker(enigma, comparator);

	if ( rank == MPI_ROOT_PROCESS_RANK ) {
		uint *messageToDecode = new uint[ MESSAGE_LENGTH ];
		generateRandomMessage(messageToDecode, MESSAGE_LENGTH);
		breaker->setMessageToDecode(messageToDecode, MESSAGE_LENGTH);

		uint *expectedMessage = new uint[ EXPECTED_MESSAGE_LENGTH ];
		copy(messageToDecode, MESSAGE_LENGTH, expectedMessage, EXPECTED_MESSAGE_LENGTH);

		breaker->setSampleToFind(expectedMessage, EXPECTED_MESSAGE_LENGTH);
	}

	cout << "Crack message - START" << endl;
	breaker->crackMessage();
	cout << "Crack message - DONE" << endl;

	if ( rank == MPI_ROOT_PROCESS_RANK ) {
		uint *result = new uint[ ROTORS ];
		breaker->getResult(result);
		cout << "Result:" << endl;
		showUint( result, ROTORS );
	}

	MPI_Finalize();

	return 0;
}

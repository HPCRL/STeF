#include "../inc/reader.h"

int main(int argc, char** argv)
{
	csf* t = (csf *) malloc(sizeof(csf));
	read_tensor(argv[1],t);
	return 0;
}
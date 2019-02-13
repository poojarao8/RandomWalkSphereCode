#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include "main.h"


// Reset the random number generator with the system clock.
void seed()
{
    srand(time(0));
}

double unifRand()
{
    return rand() / double(RAND_MAX);
}

double normalRand()
{
  std::random_device rd;
  float sample;
  // Mersenne twister PRNG, initialized with seed from previous random device instance
  std::mt19937 gen(rd());
  std::normal_distribution<float> d(0.0, 1.0);
  sample = d(gen);
  return sample;

}


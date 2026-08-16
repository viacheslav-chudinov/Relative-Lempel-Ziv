CXX = g++
CXXFLAGS = -Wall -std=c++17 -O2

TARGET = main

SRCS = $(wildcard *.cpp) ../find_a_prime/primes.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
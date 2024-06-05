# Variáveis
CXX = g++
CXXFLAGS = -Wall -I/usr/include
LDFLAGS = -L/usr/lib -lgmpxx -lgmp

# Arquivos fonte
SRCS = tp1.cpp
# Arquivos objeto
OBJS = $(SRCS:.cpp=.o)

# Nome do executável
EXEC = main

# Regra padrão
all: $(EXEC)

# Regra para criar o executável
$(EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $(EXEC) $(LDFLAGS)

# Regra para compilar os arquivos .cpp em .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Limpeza dos arquivos objetos e executáveis
clean:
	rm -f $(OBJS) $(EXEC)

# Phony targets
.PHONY: all clean

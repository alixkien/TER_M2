# Compilateur utilisé
CXX=g++

# Options communes de compilation
COMMON_FLAGS = -I ../Eigen/Eigen -std=c++11 -Wall
CXX_FLAGS = $(COMMON_FLAGS)

# Pour compiler en mode debug, lancez:
# make DEBUG=1
#####

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	# Options en mode debug - La variable est DEBUG est definie comme vraie
	CXXFLAGS += -g -Wall -DDEBUG
else
	# Options en mode optimise - La variable DEBUG est definie comme fausse
	CXXFLAGS += -O3 -std=c++11 -DNDEBUG
endif

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler
SRC = interpolation.cpp adapt.cpp ressort.cpp fonctions.cpp
main_SRC = $(SRC) main.cpp
OBJ = (main_SRC:.cpp=.o)

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : message $(main_OBJ)
	$(CXX) $(CXX_FLAGS) $(OBJ) -o $(PROG)
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) 
# et les fichiers temporaires de sauvegarde (~)
clean :
	rm -f *~ $(PROG) $(OBJ)

# Regle commune pour compiler chaque .cpp en .o
# $@ signifie "la destination" (a gauche de la regle : le fichier.o)
# $< signifie "la source" (a droite de la regle : le fichier.cpp)
%.o: %.cpp %.h
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

.PHONY: message
message:
	@echo "******************************"
	@if [ $(DEBUG) -eq "1" ]; then echo "** Building in debug mode **";  else echo "** Building in release mode **" ; fi
	@echo "******************************"


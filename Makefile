##################################################################
##################################################################
####   MAKEFILE FOR TPHIV
##################################################################
##################################################################

# ==================================
#     === Compiling Variables ===
# ==================================

C_GPP   := g++
C_CLANG := clang++

C_GPP_OPT   := -Wall -O3 -std=c++11 -Wno-predefined-identifier-outside-function -Wno-sign-compare
C_CLANG_OPT := -Wall -O3 -std=c++11 -Wno-predefined-identifier-outside-function
# -Wno-sign-compare
# -Wno-deprecated

# =================
# === Libraries ===
# =================

LINK_LIB :=        #-lgsl #-lool  # <---- For Apple Mac
LINK_LIB_GPP :=    #-lgsl -lgslcblas # <--- For Earnservs


# ==================
# === FILE LISTS ===
# ==================

SOURCE_LIST := areaUnit.cpp build_world.cpp population.cpp\
discrete_prob_dist.cpp \
simulator.cpp dcDataFrame.cpp\
dcMatrix.cpp dcTools.cpp\
disease.cpp globalvar.cpp\
intervention.cpp\
individual.cpp modelParam.cpp\
schedule.cpp socialPlace.cpp\
tests.cpp

OBJ_LIST2 := main.cpp $(SOURCE_LIST)


# ===================
# === EXECUTABLES ===
# ===================

PROG_NAME := naf
PROG_NAME_GPP := naf_gpp


# ======================
# ==== COMPILE ONLY ====
# ======================

# ==== CLANG COMPLIER ====

$(PROG_NAME): $(OBJ_LIST2)
	$(C_CLANG) $(C_CLANG_OPT) $(OBJ_LIST2) $(LINK_LIB) -o $@


# ==== GPP COMPILER (FOR EARNSERV, SHARCNET) ====

$(PROG_NAME_GPP): $(OBJ_LIST2)
	$(C_GPP) $(C_GPP_OPT) $(OBJ_LIST2) $(LINK_LIB_GPP) -o $@


# =================
# === CLEAN-UP ===
# =================

cleanout:
	rm -f *.out
	rm -f ./OUT/*.out


clean:
	rm -f *.o
	rm -f $(PROG_NAME)

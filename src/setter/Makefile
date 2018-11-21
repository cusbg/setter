# Compiler flags...
CPP_COMPILER = g++
C_COMPILER = gcc
# Include paths...
#Debug_Include_Path=-I"../../dependencies/boost" -I"../../dependencies/tbb/include" 
#Release_Include_Path=-I"../../dependencies/boost" -I"../../dependencies/tbb/include" 
#TBB include paths should be stored in the environment variable CPATH
Debug_Include_Path=-I"../../dependencies/boost" 
Release_Include_Path=-I"../../dependencies/boost" 


# Library paths...
#Debug_Library_Path=-L"../../dependencies/boost/gcclib" -L"../../dependencies/tbb/dependencies/ia32/gccvc10" 
#Release_Library_Path=-L"../../dependencies/boost/gcclib" -L"../../dependencies/tbb/dependencies/ia32/gccvc10" 
#TBB libraries paths should be in the environment variables LIBRARY_PATH and LD_LIBRARY_PATH
Debug_Library_Path=-L"../../dependencies/boost/lib" 
Release_Library_Path=-L"../../dependencies/boost/lib" 

# Additional libraries...
Debug_Libraries=-lboost_system -lboost_filesystem -lboost_date_time -lboost_program_options -lboost_serialization -Wl,--start-group -ltbb_debug  -Wl,--end-group 
Release_Libraries=-lboost_system -lboost_filesystem -lboost_date_time -lboost_program_options -lboost_serialization -Wl,--start-group -ltbb  -Wl,--end-group 

# Preprocessor definitions...
Debug_Preprocessor_Definitions=-D GCC_BUILD 
Release_Preprocessor_Definitions=-D GCC_BUILD 

# Implictly linked object files...
Debug_Implicitly_Linked_Objects=
Release_Implicitly_Linked_Objects=

# Compiler flags...
Debug_Compiler_Flags=-O0 -g -std=c++0x
Release_Compiler_Flags=-Wall -O2 -g -std=c++0x

# Builds all configurations for this project...
.PHONY: build_all_configurations
build_all_configurations: Debug Release 

# Builds the Debug configuration...
.PHONY: Debug
Debug: create_folders gccDebug/../pdb-parser/pdb_parser.o gccDebug/dist/align.o gccDebug/dist/export.o gccDebug/dist/rmsd.o gccDebug/dist/tmscore.o gccDebug/algo.o gccDebug/common.o gccDebug/cRNAStructure.o gccDebug/match_data.o gccDebug/input.o gccDebug/main.o gccDebug/make_aver_rna.o gccDebug/merge_rna.o gccDebug/multi_align.o gccDebug/multi_align_1.o gccDebug/multi_align_2.o gccDebug/multi_align_0.o gccDebug/parallel.o gccDebug/parallel_make_aver.o gccDebug/parameters.o gccDebug/print_to_pdb.o gccDebug/rotation.o gccDebug/serializer.o 
	g++ gccDebug/../../pdb-parser/pdb_parser.o gccDebug/dist/align.o gccDebug/dist/export.o gccDebug/dist/rmsd.o gccDebug/dist/tmscore.o gccDebug/algo.o gccDebug/common.o gccDebug/cRNAStructure.o gccDebug/match_data.o gccDebug/input.o gccDebug/main.o gccDebug/make_aver_rna.o gccDebug/merge_rna.o gccDebug/multi_align.o gccDebug/multi_align_1.o gccDebug/multi_align_2.o gccDebug/multi_align_0.o gccDebug/parallel.o gccDebug/parallel_make_aver.o gccDebug/parameters.o gccDebug/print_to_pdb.o gccDebug/rotation.o gccDebug/serializer.o  $(Debug_Library_Path) $(Debug_Libraries) -Wl,-rpath,./ -o ../debug/multisetter

# Compiles file ../pdb-parser/pdb_parser.cpp for the Debug configuration...
-include gccDebug/../../pdb-parser/pdb_parser.d
gccDebug/../pdb-parser/pdb_parser.o: ../pdb-parser/pdb_parser.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c ../pdb-parser/pdb_parser.cpp $(Debug_Include_Path) -o gccDebug/../../pdb-parser/pdb_parser.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM ../pdb-parser/pdb_parser.cpp $(Debug_Include_Path) > gccDebug/../../pdb-parser/pdb_parser.d

# Compiles file dist/align.cpp for the Debug configuration...
-include gccDebug/dist/align.d
gccDebug/dist/align.o: dist/align.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c dist/align.cpp $(Debug_Include_Path) -o gccDebug/dist/align.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM dist/align.cpp $(Debug_Include_Path) > gccDebug/dist/align.d

# Compiles file dist/export.cpp for the Debug configuration...
-include gccDebug/dist/export.d
gccDebug/dist/export.o: dist/export.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c dist/export.cpp $(Debug_Include_Path) -o gccDebug/dist/export.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM dist/export.cpp $(Debug_Include_Path) > gccDebug/dist/export.d

# Compiles file dist/rmsd.cpp for the Debug configuration...
-include gccDebug/dist/rmsd.d
gccDebug/dist/rmsd.o: dist/rmsd.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c dist/rmsd.cpp $(Debug_Include_Path) -o gccDebug/dist/rmsd.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM dist/rmsd.cpp $(Debug_Include_Path) > gccDebug/dist/rmsd.d

# Compiles file dist/tmscore.cpp for the Debug configuration...
-include gccDebug/dist/tmscore.d
gccDebug/dist/tmscore.o: dist/tmscore.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c dist/tmscore.cpp $(Debug_Include_Path) -o gccDebug/dist/tmscore.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM dist/tmscore.cpp $(Debug_Include_Path) > gccDebug/dist/tmscore.d

# Compiles file algo.cpp for the Debug configuration...
-include gccDebug/algo.d
gccDebug/algo.o: algo.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c algo.cpp $(Debug_Include_Path) -o gccDebug/algo.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM algo.cpp $(Debug_Include_Path) > gccDebug/algo.d

# Compiles file common.cpp for the Debug configuration...
-include gccDebug/common.d
gccDebug/common.o: common.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c common.cpp $(Debug_Include_Path) -o gccDebug/common.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM common.cpp $(Debug_Include_Path) > gccDebug/common.d

# Compiles file cRNAStructure.cpp for the Debug configuration...
-include gccDebug/cRNAStructure.d
gccDebug/cRNAStructure.o: cRNAStructure.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c cRNAStructure.cpp $(Debug_Include_Path) -o gccDebug/cRNAStructure.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM cRNAStructure.cpp $(Debug_Include_Path) > gccDebug/cRNAStructure.d

# Compiles file match_data.cpp for the Debug configuration...
-include gccDebug/match_data.d
gccDebug/match_data.o: match_data.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c match_data.cpp $(Debug_Include_Path) -o gccDebug/match_data.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM match_data.cpp $(Debug_Include_Path) > gccDebug/match_data.d

# Compiles file input.cpp for the Debug configuration...
-include gccDebug/input.d
gccDebug/input.o: input.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c input.cpp $(Debug_Include_Path) -o gccDebug/input.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM input.cpp $(Debug_Include_Path) > gccDebug/input.d

# Compiles file main.cpp for the Debug configuration...
-include gccDebug/main.d
gccDebug/main.o: main.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c main.cpp $(Debug_Include_Path) -o gccDebug/main.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM main.cpp $(Debug_Include_Path) > gccDebug/main.d

# Compiles file make_aver_rna.cpp for the Debug configuration...
-include gccDebug/make_aver_rna.d
gccDebug/make_aver_rna.o: make_aver_rna.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c make_aver_rna.cpp $(Debug_Include_Path) -o gccDebug/make_aver_rna.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM make_aver_rna.cpp $(Debug_Include_Path) > gccDebug/make_aver_rna.d

# Compiles file merge_rna.cpp for the Debug configuration...
-include gccDebug/merge_rna.d
gccDebug/merge_rna.o: merge_rna.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c merge_rna.cpp $(Debug_Include_Path) -o gccDebug/merge_rna.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM merge_rna.cpp $(Debug_Include_Path) > gccDebug/merge_rna.d

# Compiles file multi_align.cpp for the Debug configuration...
-include gccDebug/multi_align.d
gccDebug/multi_align.o: multi_align.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c multi_align.cpp $(Debug_Include_Path) -o gccDebug/multi_align.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM multi_align.cpp $(Debug_Include_Path) > gccDebug/multi_align.d

# Compiles file multi_align_1.cpp for the Debug configuration...
-include gccDebug/multi_align_1.d
gccDebug/multi_align_1.o: multi_align_1.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c multi_align_1.cpp $(Debug_Include_Path) -o gccDebug/multi_align_1.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM multi_align_1.cpp $(Debug_Include_Path) > gccDebug/multi_align_1.d

# Compiles file multi_align_2.cpp for the Debug configuration...
-include gccDebug/multi_align_2.d
gccDebug/multi_align_2.o: multi_align_2.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c multi_align_2.cpp $(Debug_Include_Path) -o gccDebug/multi_align_2.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM multi_align_2.cpp $(Debug_Include_Path) > gccDebug/multi_align_2.d

# Compiles file multi_align_0.cpp for the Debug configuration...
-include gccDebug/multi_align_0.d
gccDebug/multi_align_0.o: multi_align_0.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c multi_align_0.cpp $(Debug_Include_Path) -o gccDebug/multi_align_0.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM multi_align_0.cpp $(Debug_Include_Path) > gccDebug/multi_align_0.d

# Compiles file parallel.cpp for the Debug configuration...
-include gccDebug/parallel.d
gccDebug/parallel.o: parallel.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c parallel.cpp $(Debug_Include_Path) -o gccDebug/parallel.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM parallel.cpp $(Debug_Include_Path) > gccDebug/parallel.d

# Compiles file parallel_make_aver.cpp for the Debug configuration...
-include gccDebug/parallel_make_aver.d
gccDebug/parallel_make_aver.o: parallel_make_aver.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c parallel_make_aver.cpp $(Debug_Include_Path) -o gccDebug/parallel_make_aver.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM parallel_make_aver.cpp $(Debug_Include_Path) > gccDebug/parallel_make_aver.d

# Compiles file parameters.cpp for the Debug configuration...
-include gccDebug/parameters.d
gccDebug/parameters.o: parameters.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c parameters.cpp $(Debug_Include_Path) -o gccDebug/parameters.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM parameters.cpp $(Debug_Include_Path) > gccDebug/parameters.d

# Compiles file print_to_pdb.cpp for the Debug configuration...
-include gccDebug/print_to_pdb.d
gccDebug/print_to_pdb.o: print_to_pdb.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c print_to_pdb.cpp $(Debug_Include_Path) -o gccDebug/print_to_pdb.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM print_to_pdb.cpp $(Debug_Include_Path) > gccDebug/print_to_pdb.d

# Compiles file rna_structure_serialization.cpp for the Debug configuration...
-include gccDebug/rna_structure_serialization.d
gccDebug/rna_structure_serialization.o: rna_structure_serialization.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c rna_structure_serialization.cpp $(Debug_Include_Path) -o gccDebug/rna_structure_serialization.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM rna_structure_serialization.cpp $(Debug_Include_Path) > gccDebug/rna_structure_serialization.d

# Compiles file rotation.cpp for the Debug configuration...
-include gccDebug/rotation.d
gccDebug/rotation.o: rotation.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c rotation.cpp $(Debug_Include_Path) -o gccDebug/rotation.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM rotation.cpp $(Debug_Include_Path) > gccDebug/rotation.d

# Compiles file serializer.cpp for the Debug configuration...
-include gccDebug/serializer.d
gccDebug/serializer.o: serializer.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c serializer.cpp $(Debug_Include_Path) -o gccDebug/serializer.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM serializer.cpp $(Debug_Include_Path) > gccDebug/serializer.d

# Builds the Release configuration...
.PHONY: Release
Release: create_folders gccRelease/../../pdb-parser/pdb_parser.o gccRelease/dist/align.o gccRelease/dist/export.o gccRelease/dist/rmsd.o gccRelease/dist/tmscore.o gccRelease/algo.o gccRelease/common.o gccRelease/cRNAStructure.o gccRelease/match_data.o gccRelease/input.o gccRelease/main.o gccRelease/make_aver_rna.o gccRelease/merge_rna.o gccRelease/multi_align.o gccRelease/multi_align_1.o gccRelease/multi_align_2.o gccRelease/multi_align_0.o gccRelease/parallel.o gccRelease/parallel_make_aver.o gccRelease/parameters.o gccRelease/print_to_pdb.o gccRelease/rotation.o gccRelease/serializer.o 
	g++ gccRelease/../../pdb-parser/pdb_parser.o gccRelease/dist/align.o gccRelease/dist/export.o gccRelease/dist/rmsd.o gccRelease/dist/tmscore.o gccRelease/algo.o gccRelease/common.o gccRelease/cRNAStructure.o gccRelease/match_data.o gccRelease/input.o gccRelease/main.o gccRelease/make_aver_rna.o gccRelease/merge_rna.o gccRelease/multi_align.o gccRelease/multi_align_1.o gccRelease/multi_align_2.o gccRelease/multi_align_0.o gccRelease/parallel.o gccRelease/parallel_make_aver.o gccRelease/parameters.o gccRelease/print_to_pdb.o gccRelease/rotation.o gccRelease/serializer.o  $(Release_Library_Path) $(Release_Libraries) -Wl,-rpath,./ -o ../release/multisetter

# Compiles file ../pdb-parser/pdb_parser.cpp for the Release configuration...
-include gccRelease/../../pdb-parser/pdb_parser.d
gccRelease/../../pdb-parser/pdb_parser.o: ../pdb-parser/pdb_parser.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c ../pdb-parser/pdb_parser.cpp $(Release_Include_Path) -o gccRelease/../../pdb-parser/pdb_parser.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM ../pdb-parser/pdb_parser.cpp $(Release_Include_Path) > gccRelease/../../pdb-parser/pdb_parser.d

# Compiles file dist/align.cpp for the Release configuration...
-include gccRelease/dist/align.d
gccRelease/dist/align.o: dist/align.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c dist/align.cpp $(Release_Include_Path) -o gccRelease/dist/align.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM dist/align.cpp $(Release_Include_Path) > gccRelease/dist/align.d

# Compiles file dist/export.cpp for the Release configuration...
-include gccRelease/dist/export.d
gccRelease/dist/export.o: dist/export.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c dist/export.cpp $(Release_Include_Path) -o gccRelease/dist/export.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM dist/export.cpp $(Release_Include_Path) > gccRelease/dist/export.d

# Compiles file dist/rmsd.cpp for the Release configuration...
-include gccRelease/dist/rmsd.d
gccRelease/dist/rmsd.o: dist/rmsd.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c dist/rmsd.cpp $(Release_Include_Path) -o gccRelease/dist/rmsd.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM dist/rmsd.cpp $(Release_Include_Path) > gccRelease/dist/rmsd.d

# Compiles file dist/tmscore.cpp for the Release configuration...
-include gccRelease/dist/tmscore.d
gccRelease/dist/tmscore.o: dist/tmscore.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c dist/tmscore.cpp $(Release_Include_Path) -o gccRelease/dist/tmscore.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM dist/tmscore.cpp $(Release_Include_Path) > gccRelease/dist/tmscore.d

# Compiles file algo.cpp for the Release configuration...
-include gccRelease/algo.d
gccRelease/algo.o: algo.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c algo.cpp $(Release_Include_Path) -o gccRelease/algo.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM algo.cpp $(Release_Include_Path) > gccRelease/algo.d

# Compiles file common.cpp for the Release configuration...
-include gccRelease/common.d
gccRelease/common.o: common.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c common.cpp $(Release_Include_Path) -o gccRelease/common.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM common.cpp $(Release_Include_Path) > gccRelease/common.d

# Compiles file cRNAStructure.cpp for the Release configuration...
-include gccRelease/cRNAStructure.d
gccRelease/cRNAStructure.o: cRNAStructure.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c cRNAStructure.cpp $(Release_Include_Path) -o gccRelease/cRNAStructure.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM cRNAStructure.cpp $(Release_Include_Path) > gccRelease/cRNAStructure.d

# Compiles file match_data.cpp for the Release configuration...
-include gccRelease/match_data.d
gccRelease/match_data.o: match_data.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c match_data.cpp $(Release_Include_Path) -o gccRelease/match_data.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM match_data.cpp $(Release_Include_Path) > gccRelease/match_data.d

# Compiles file input.cpp for the Release configuration...
-include gccRelease/input.d
gccRelease/input.o: input.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c input.cpp $(Release_Include_Path) -o gccRelease/input.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM input.cpp $(Release_Include_Path) > gccRelease/input.d

# Compiles file main.cpp for the Release configuration...
-include gccRelease/main.d
gccRelease/main.o: main.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c main.cpp $(Release_Include_Path) -o gccRelease/main.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM main.cpp $(Release_Include_Path) > gccRelease/main.d

# Compiles file make_aver_rna.cpp for the Release configuration...
-include gccRelease/make_aver_rna.d
gccRelease/make_aver_rna.o: make_aver_rna.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c make_aver_rna.cpp $(Release_Include_Path) -o gccRelease/make_aver_rna.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM make_aver_rna.cpp $(Release_Include_Path) > gccRelease/make_aver_rna.d

# Compiles file merge_rna.cpp for the Release configuration...
-include gccRelease/merge_rna.d
gccRelease/merge_rna.o: merge_rna.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c merge_rna.cpp $(Release_Include_Path) -o gccRelease/merge_rna.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM merge_rna.cpp $(Release_Include_Path) > gccRelease/merge_rna.d

# Compiles file multi_align.cpp for the Release configuration...
-include gccRelease/multi_align.d
gccRelease/multi_align.o: multi_align.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c multi_align.cpp $(Release_Include_Path) -o gccRelease/multi_align.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM multi_align.cpp $(Release_Include_Path) > gccRelease/multi_align.d

# Compiles file multi_align_1.cpp for the Release configuration...
-include gccRelease/multi_align_1.d
gccRelease/multi_align_1.o: multi_align_1.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c multi_align_1.cpp $(Release_Include_Path) -o gccRelease/multi_align_1.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM multi_align_1.cpp $(Release_Include_Path) > gccRelease/multi_align_1.d

# Compiles file multi_align_2.cpp for the Release configuration...
-include gccRelease/multi_align_2.d
gccRelease/multi_align_2.o: multi_align_2.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c multi_align_2.cpp $(Release_Include_Path) -o gccRelease/multi_align_2.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM multi_align_2.cpp $(Release_Include_Path) > gccRelease/multi_align_2.d

# Compiles file multi_align_0.cpp for the Release configuration...
-include gccRelease/multi_align_0.d
gccRelease/multi_align_0.o: multi_align_0.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c multi_align_0.cpp $(Release_Include_Path) -o gccRelease/multi_align_0.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM multi_align_0.cpp $(Release_Include_Path) > gccRelease/multi_align_0.d

# Compiles file parallel.cpp for the Release configuration...
-include gccRelease/parallel.d
gccRelease/parallel.o: parallel.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c parallel.cpp $(Release_Include_Path) -o gccRelease/parallel.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM parallel.cpp $(Release_Include_Path) > gccRelease/parallel.d

# Compiles file parallel_make_aver.cpp for the Release configuration...
-include gccRelease/parallel_make_aver.d
gccRelease/parallel_make_aver.o: parallel_make_aver.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c parallel_make_aver.cpp $(Release_Include_Path) -o gccRelease/parallel_make_aver.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM parallel_make_aver.cpp $(Release_Include_Path) > gccRelease/parallel_make_aver.d

# Compiles file parameters.cpp for the Release configuration...
-include gccRelease/parameters.d
gccRelease/parameters.o: parameters.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c parameters.cpp $(Release_Include_Path) -o gccRelease/parameters.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM parameters.cpp $(Release_Include_Path) > gccRelease/parameters.d

# Compiles file print_to_pdb.cpp for the Release configuration...
-include gccRelease/print_to_pdb.d
gccRelease/print_to_pdb.o: print_to_pdb.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c print_to_pdb.cpp $(Release_Include_Path) -o gccRelease/print_to_pdb.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM print_to_pdb.cpp $(Release_Include_Path) > gccRelease/print_to_pdb.d

# Compiles file rna_structure_serialization.cpp for the Release configuration...
-include gccRelease/rna_structure_serialization.d
gccRelease/rna_structure_serialization.o: rna_structure_serialization.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c rna_structure_serialization.cpp $(Release_Include_Path) -o gccRelease/rna_structure_serialization.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM rna_structure_serialization.cpp $(Release_Include_Path) > gccRelease/rna_structure_serialization.d

# Compiles file rotation.cpp for the Release configuration...
-include gccRelease/rotation.d
gccRelease/rotation.o: rotation.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c rotation.cpp $(Release_Include_Path) -o gccRelease/rotation.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM rotation.cpp $(Release_Include_Path) > gccRelease/rotation.d

# Compiles file serializer.cpp for the Release configuration...
-include gccRelease/serializer.d
gccRelease/serializer.o: serializer.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c serializer.cpp $(Release_Include_Path) -o gccRelease/serializer.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM serializer.cpp $(Release_Include_Path) > gccRelease/serializer.d

# Creates the intermediate and output folders for each configuration...
.PHONY: create_folders
create_folders:
	mkdir -p gccDebug/source
	mkdir -p gccDebug/dist
	mkdir -p ../debug
	mkdir -p gccRelease/source
	mkdir -p gccRelease/dist
	mkdir -p ../release



# Cleans intermediate and output files (objects, libraries, executables)...
.PHONY: clean
clean:
	rm -f gccDebug/*.o
	rm -f gccDebug/*.d
	rm -f ../debug/*
	rm -f gccRelease/*.o
	rm -f gccRelease/*.d
	rm -f ../release/*


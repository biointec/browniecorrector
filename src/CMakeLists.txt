add_executable(brownie graphaln.cpp nodechain.cpp kmercounttable.cpp kmeroverlaptable.cpp readcorrection.cpp alignment.cpp coverage.cpp library.cpp kmernpp.cpp kmertable.cpp correctgraph.cpp dsnode.cpp nucleotide.cpp nodeendstable.cpp settings.cpp util.cpp tstring.cpp kmeroverlap.cpp graph.cpp brownie.cpp refcomp.cpp suffix_tree.c tipClipping.cpp markovChain.cpp bubbleDetection.cpp flowCorrection.cpp)
add_executable(louvain louvainSrc/main_louvain.cpp louvainSrc/zahn.cpp louvainSrc/balmod.cpp louvainSrc/condora.cpp louvainSrc/devind.cpp  louvainSrc/devuni.cpp louvainSrc/dp.cpp louvainSrc/goldberg.cpp louvainSrc/graph_binary.cpp  louvainSrc/louvain.cpp  louvainSrc/modularity.cpp louvainSrc/owzad.cpp louvainSrc/quality.cpp louvainSrc/shimalik.cpp )

add_executable(readclustering readClusteringSrc/main.cpp readClusteringSrc/connectedComponent.cpp readClusteringSrc/calSim.cpp util.cpp readClusteringSrc/extractReads.cpp library.cpp nucleotide.cpp readClusteringSrc/kmerObserver.cpp )

add_executable(hierarchy louvainSrc/main_hierarchy.cpp   )

add_executable(convert louvainSrc/main_convert.cpp  louvainSrc/graph.cpp )

add_executable(matrix louvainSrc/main_matrix.cpp  )

#add_executable(calSim readClusteringSrc/all2allOL.cpp)

target_link_libraries(brownie readfile essaMEM pthread)
target_link_libraries(readclustering readfile)
target_link_libraries(readclustering pthread)

if (ZLIB_FOUND)
   target_link_libraries(brownie ${ZLIB_LIBRARY})
   target_link_libraries(readclustering ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)

install(TARGETS brownie RUNTIME DESTINATION bin)

add_subdirectory(readfile)
add_subdirectory(essaMEM-master)



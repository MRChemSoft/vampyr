target_sources(vampyr3d
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/operators.cpp
	${CMAKE_CURRENT_LIST_DIR}/ABGVOperator.cpp
	${CMAKE_CURRENT_LIST_DIR}/BSOperator.cpp
	${CMAKE_CURRENT_LIST_DIR}/ConvolutionOperator.cpp
	${CMAKE_CURRENT_LIST_DIR}/DerivativeKernel.cpp
	${CMAKE_CURRENT_LIST_DIR}/DerivativeOperator.cpp
	${CMAKE_CURRENT_LIST_DIR}/HelmholtzKernel.cpp
	${CMAKE_CURRENT_LIST_DIR}/HelmholtzOperator.cpp
	${CMAKE_CURRENT_LIST_DIR}/IdentityConvolution.cpp
	${CMAKE_CURRENT_LIST_DIR}/IdentityKernel.cpp
	${CMAKE_CURRENT_LIST_DIR}/MWOperator.cpp
	${CMAKE_CURRENT_LIST_DIR}/OperatorState.cpp
	${CMAKE_CURRENT_LIST_DIR}/OperatorStatistics.cpp
	${CMAKE_CURRENT_LIST_DIR}/PHOperator.cpp
	${CMAKE_CURRENT_LIST_DIR}/PoissonOperator.cpp

  )





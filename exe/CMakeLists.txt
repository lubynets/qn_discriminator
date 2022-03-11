# add_executable(discriminator discriminator.cxx)
# add_dependencies(discriminator QnAnalysisDiscriminator)
# target_link_libraries(discriminator QnAnalysisDiscriminator)
# 
# add_executable(shapewriter shapewriter.cxx)
# add_dependencies(shapewriter QnAnalysisDiscriminator)
# target_link_libraries(shapewriter QnAnalysisDiscriminator)
# 
# add_executable(shapetensorwriter shapetensorwriter.cxx)
# add_dependencies(shapetensorwriter QnAnalysisDiscriminator)
# target_link_libraries(shapetensorwriter QnAnalysisDiscriminator)

add_executable(shapetensorwriter_fit shapetensorwriter_fit.cxx)
add_dependencies(shapetensorwriter_fit QnDiscriminator)
target_link_libraries(shapetensorwriter_fit QnTools::Base ${ROOT_LIBRARIES} QnDiscriminator)

# add_executable(discriminatortensor discriminatortensor.cxx)
# add_dependencies(discriminatortensor QnAnalysisDiscriminator)
# target_link_libraries(discriminatortensor QnAnalysisDiscriminator)
# 
# add_executable(mcfitter mcfitter.cxx)
# add_dependencies(mcfitter QnAnalysisDiscriminator)
# target_link_libraries(mcfitter QnAnalysisDiscriminator)
# 
# add_executable(mcextractor mcextractor.cxx)
# add_dependencies(mcextractor QnAnalysisDiscriminator)
# target_link_libraries(mcextractor QnAnalysisDiscriminator)
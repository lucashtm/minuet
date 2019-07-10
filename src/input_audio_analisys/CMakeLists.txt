
add_executable(input_audio_analisys AudioFile.cpp signal_analyser.cpp)
target_link_libraries(input_audio_analisys Qt5::Core Qt5::Gui Qt5::Multimedia Qt5::Widgets -lfftw3) 
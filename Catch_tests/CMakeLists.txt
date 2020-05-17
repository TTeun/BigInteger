add_executable(Catch_tests_run BigUIntTest.cpp ../BigInt/BigUInt.cpp ../BigInt/BigUInt.h ../timer.cpp ../timer.h ../BigInt/DigitVector.cpp ../BigInt/DigitVector.h ../BigInt/BigUIntBase.cpp ../BigInt/BigUIntBase.h ../BigInt/BigInt.cpp ../BigInt/BigInt.h)

target_link_libraries(Catch_tests_run)


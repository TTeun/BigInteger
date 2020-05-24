#define CATCH_CONFIG_MAIN

#include "../BigInt/BigUInt.h"

#include "catch.hpp"

using namespace big;

TEST_CASE("BigUInt tests", "[BigUIntTests]") {

    SECTION("Sum of natural numbers") {
        BigUInt sumOfNaturalNumbers("5000050000");
        BigUInt b(0);
        for (size_t i = 0; i != 100001; ++i) { b += BigUInt(i); }
        CHECK(b == sumOfNaturalNumbers);
        for (size_t i = 0; i != 100001; ++i) { b -= BigUInt(i); }
        CHECK(b == 0);
    }

    SECTION("Modulo") {
        BigUInt b("9945549994547342487362847632487678463291881729");
        CHECK(b % 1237862343UL == 688386621UL);
    }

    SECTION("Fixed division") {
        const auto a = BigUInt("31988428769874667983746348976357643987568732476328746237847832647832676327573"
                               "28657832647832647832784"
                               "683724678326487326743656329847047017259887150872142");
        const auto b = BigUInt("29371982479821749842772102198749275124736283768732687126321");
        const auto c = BigUInt("108907966262918718058705949100715069245292669239739335359987042442051400703079828741589494778");
        CHECK(a / b == c);
    }
    SECTION("Fixed subtraction") {
        BigUInt b("21784672876238762183782176387213213");
        BigUInt c("327846728238762183782176387213213");
        b -= c;
        assert(b == BigUInt("21456826148000000000000000000000000"));
    }

    SECTION("Addition") {
        for (size_t dummy = 0; dummy != 200; ++dummy) {
            size_t a = rand() % (std::numeric_limits<size_t>::max() / 2ul);
            size_t b = rand() % (std::numeric_limits<size_t>::max() / 2ul);
            const auto k = BigUInt(a) + BigUInt(b);
            CHECK(k == a + b);
        }
    }

    SECTION("Multiplication") {
        for (size_t dummy = 0; dummy != 200; ++dummy) {
            size_t a = rand() % std::numeric_limits<size_t>::max();
            size_t b = rand() % std::numeric_limits<size_t>::max();
            const auto k = BigUInt(a) * BigUInt(b);
            CHECK(BigUInt(a) * BigUInt(b) == a * b);
        }
    }

    SECTION("Division") {
        for (size_t dummy = 0; dummy != 200; ++dummy) {
            const auto a = static_cast<size_t>(rand() + 1ul);
            const auto b = static_cast<size_t>(rand() + 1ul);
            const auto k = BigUInt(a) / BigUInt(b);
            CHECK(k == a / b);
        }
    }

    SECTION("Power") {
        const BigUInt twoToThePowerTwoThousand("114813069527425452423283320117768198402231770208869520047764273682"
                                               "576626139237031385665948631650626991"
                                               "844596463898746277344711896086305533142593135616665318539129989145"
                                               "312280000688779148240044871428926990"
                                               "063486244781615463646388363947317026040466353970904996558162398808"
                                               "944629605623311649536164221970332681"
                                               "344168908984458505602379484807914058900934776500429002716706625830"
                                               "522008132236281291761267883317206598"
                                               "995396418127021779858404042159853183251540889433902091920554957783"
                                               "589672039160081957216630582755380425"
                                               "583726015528348786419432054508915275783882625175435528800822842770"
                                               "817965453762184851149029376");
        BigUInt b = power(BigUInt(2), 2000);
        CHECK(b == twoToThePowerTwoThousand);
    }

    SECTION("Factorial") {
        const BigUInt thousandFactorial("402387260077093773543702433923003985719374864210714632543799910429"
                                        "938512398629020592044208486969404"
                                        "800479988610197196058631666872994808558901323829669944590997424504"
                                        "087073759918823627727188732519779"
                                        "505950995276120874975462497043601418278094646496291056393887437886"
                                        "487337119181045825783647849977012"
                                        "476632889835955735432513185323958463075557409114262417474349347553"
                                        "428646576611667797396668820291207"
                                        "379143853719588249808126867838374559731746136085379534524221586593"
                                        "201928090878297308431392844403281"
                                        "231558611036976801357304216168747609675871348312025478589320767169"
                                        "132448426236131412508780208000261"
                                        "683151027341827977704784635868170164365024153691398281264810213092"
                                        "761244896359928705114964975419909"
                                        "342221566832572080821333186116811553615836546984046708975602900950"
                                        "537616475847728421889679646244945"
                                        "160765353408198901385442487984959953319101723355556602139450399736"
                                        "280750137837615307127761926849034"
                                        "352625200015888535147331611702103968175921510907788019393178114194"
                                        "545257223865541461062892187960223"
                                        "838971476088506276862967146674697562911234082439208160153780889893"
                                        "964518263243671616762179168909779"
                                        "911903754031274622289988005195444414282012187361745992642956581746"
                                        "628302955570299024324153181617210"
                                        "465832036786906117260158783520751516284225540265170483304226143974"
                                        "286933061690897968482590125458327"
                                        "168226458066526769958652682272807075781391858178889652208164348344"
                                        "825993266043367660176999612831860"
                                        "788386150279465955131156552036093988180612138558600301435694527224"
                                        "206344631797460594682573103790084"
                                        "024432438465657245014402821885252470935190620929023136493273497565"
                                        "513958720559654228749774011413346"
                                        "962715422845862377387538230483865688976461927383814900140767310446"
                                        "640259899490222221765904339901886"
                                        "018566526485061799702356193897017860040811889729918311021171229845"
                                        "901641921068884387121855646124960"
                                        "798722908519296819372388642614839657382291123125024186649353143970"
                                        "137428531926649875337218940694281"
                                        "434118520158014123344828015051399694290153483077644569099073152433"
                                        "278288269864602789864321139083506"
                                        "217095002597389863554277196742822248757586765752344220207573630569"
                                        "498825087968928162753848863396909"
                                        "959826280956121450994871701244516461260379029309120889086942028510"
                                        "640182154399457156805941872748998"
                                        "094254742173582401063677404595741785160829230135358081840096996372"
                                        "524230560855903700624271243416909"
                                        "004153690105933983835777939410970027753472000000000000000000000000"
                                        "000000000000000000000000000000000"
                                        "000000000000000000000000000000000000000000000000000000000000000000"
                                        "000000000000000000000000000000000"
                                        "000000000000000000000000000000000000000000000000000000000000000000"
                                        "000000000000000000000000000");
        BigUInt b = 1;
        for (size_t i = 1; i != 1001; ++i) { b *= i; }
        CHECK(b == thousandFactorial);
    }
}
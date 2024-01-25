#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>

#include "../Helpers/Assert.h"
#include "../Helpers/Debug.h"
#include "../Helpers/Timer.h"
#include "../Helpers/Vector/Vector.h"

#include "../DataStructures/Graph/Classes/DynamicGraph.h"
#include "../DataStructures/Graph/Classes/StaticGraph.h"
#include "../DataStructures/Graph/Classes/EdgeList.h"

#include "../DataStructures/Graph/Graph.h"

#include "../UnitTests/PublicTransitProfiles/CrossTransfer.h"
#include "../UnitTests/PublicTransitProfiles/JourneyOrder.h"
#include "../UnitTests/PublicTransitProfiles/TransferLoops.h"
#include "../UnitTests/PublicTransitProfiles/UndominatedWalking.h"
#include "../UnitTests/PublicTransitProfiles/UnrestrictedWalking.h"

#include "../UnitTests/Partition/NestedDissection.h"
#include "../UnitTests/Partition/Dinic.h"
#include "../UnitTests/Partition/MinimumBipartiteVertexCover.h"

#include "../UnitTests/RAPTOR/ReverseNetwork.h"

#include "../UnitTests/Graph/DynamicGraph.h"
#include "../UnitTests/Graph/StaticGraph.h"
#include "../UnitTests/Graph/EdgeList.h"
#include "../UnitTests/Graph/UnionAndIntersection.h"

#include "../UnitTests/Attributes/Attributes.h"

#include "../UnitTests/DataStructures/CoordinateTree.h"

#include "../UnitTests/Helpers/IO.h"
#include "../UnitTests/Helpers/Permutation.h"
#include "../UnitTests/Helpers/FileSystem.h"

#include "../Algorithms/CSA/OldProfileCSA.h"
#include "../Algorithms/RAPTOR/AlternatingRAPTOR.h"
#include "../Algorithms/RAPTOR/RangeRAPTOR/RangeRAPTOR.h"
#include "../UnitTests/RouteFlags/RouteFlagsPreprocessing.h"
#include "../UnitTests/RouteFlags/RouteFlagsQuery.h"

#include "../UnitTests/ULTRA/ULTRATests.h"

using RangeRAPTOR = RAPTOR::RangeRAPTOR::RangeDijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, false, true, true, true, true>;
using AlternatingRAPTOR = RAPTOR::AlternatingDijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, false, true>;
using ProfileCSA = CSA::OldProfileCSA<12, true, true, false>;

int main(int, char**) {
    checkAsserts();

    UnitTests::Attributes().check();

    UnitTests::CoordinateTree().check();

    UnitTests::IO().check();
    UnitTests::Permutation().check();
    UnitTests::FileSystem().check();

    UnitTests::DynamicGraph().check();
    UnitTests::StaticGraph().check();
    UnitTests::EdgeList().check();
    UnitTests::UnionAndIntersection().check();

    UnitTests::TransferLoops().checkRAPTOR<RangeRAPTOR>("RangeRAPTOR");
    UnitTests::TransferLoops().checkRAPTOR<AlternatingRAPTOR>("AlternatingRAPTOR");
    UnitTests::TransferLoops().checkCSA<ProfileCSA>("ProfileCSA");

    UnitTests::CrossTransfers().checkRAPTOR<RangeRAPTOR>("RangeRAPTOR");
    UnitTests::CrossTransfers().checkRAPTOR<AlternatingRAPTOR>("AlternatingRAPTOR");
    UnitTests::CrossTransfers().checkCSA<ProfileCSA>("ProfileCSA");

    UnitTests::UndominatedWalking().checkRAPTOR<RangeRAPTOR>("RangeRAPTOR");
    UnitTests::UndominatedWalking().checkRAPTOR<AlternatingRAPTOR>("AlternatingRAPTOR");
    UnitTests::UndominatedWalking().checkCSA<ProfileCSA>("ProfileCSA");

    UnitTests::JourneyOrder().checkRAPTOR<RangeRAPTOR>("RangeRAPTOR");
    UnitTests::JourneyOrder().checkRAPTOR<AlternatingRAPTOR>("AlternatingRAPTOR");
    UnitTests::JourneyOrder().checkCSA<ProfileCSA>("ProfileCSA");

    UnitTests::UnrestrictedWalking().checkRAPTOR<RangeRAPTOR>("RangeRAPTOR");
    UnitTests::UnrestrictedWalking().checkRAPTOR<AlternatingRAPTOR>("AlternatingRAPTOR");
    UnitTests::UnrestrictedWalking().checkCSA<ProfileCSA>("ProfileCSA");

    UnitTests::ReverseNetwork().check<RangeRAPTOR>();

    UnitTests::NestedDissection().check();
    UnitTests::Dinic().check();
    UnitTests::MinimumBipartiteVertexCover().check();

    UnitTests::RouteFlagsPreprocessing().check();

    UnitTests::RouteFlagsQuery<RAPTOR::RouteFlags::StopEventFlags>().check();
    UnitTests::RouteFlagsQuery<RAPTOR::RouteFlags::TripFlags>().check();
    UnitTests::RouteFlagsQuery<RAPTOR::RouteFlags::RouteFlags>().check();

    UnitTests::ULTRATests().checkWeakDominationStopToStop();
    UnitTests::ULTRATests().checkWeakDominationEventToEvent();
    UnitTests::ULTRATests().checkCyclicalDominationStopToStop();
    UnitTests::ULTRATests().checkCyclicalDominationEventToEvent();
    UnitTests::ULTRATests().checkRouteOrderStopToStop();
    UnitTests::ULTRATests().checkRouteOrderEventToEvent();

    UnitTests::evaluate();

    return 0;
}

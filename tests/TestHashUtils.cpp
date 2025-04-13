#define CATCH_CONFIG_MAIN

#include "../include/HashUtils.h"
#include "../include/Types.h"
#include "../include/external/MurmurHash3.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE("initPairwiseHashes generates valid hash seeds", "[HashUtils]") {
  initPairwiseHashes();

  auto &ctx = HashContext::instance();
  REQUIRE(ctx.seed1 >= 1);
  REQUIRE(ctx.seed2 >= 1);
  REQUIRE(ctx.seed1 <= PRIME - 1);
  REQUIRE(ctx.seed2 <= PRIME - 1);
  REQUIRE(ctx.seed1 != ctx.seed2);
}

TEST_CASE("hashAC combines h1(a) and h2(c) correctly", "[HashUtils]") {
  double comb1 = hashAC(0.0, 0.0);
  double comb2 = hashAC(0.1, 0.9);
  double comb3 = hashAC(0.99999, 0.00001);

  REQUIRE(comb1 >= 0.0);
  REQUIRE(comb1 <= 1.0);
  REQUIRE(comb2 >= 0.0);
  REQUIRE(comb2 <= 1.0);
  REQUIRE(comb3 >= 0.0);
  REQUIRE(comb3 <= 1.0);

  // Check for variety and that results differ with inputs [0,1)
  REQUIRE(hashAC(0.0, 0.7) != hashAC(0.0, 0.0));
  REQUIRE(hashAC(0.5, 0.5) == hashAC(0.4, 0.4));
  REQUIRE(hashAC(0.5, 0.5) == hashAC(0.5, 0.5));
  REQUIRE(hashAC(0.6, 0.5) != hashAC(0.5, 0.6));
}

TEST_CASE("murmur_hash returns values in [0, 1]", "[HashUtils]") {
  initPairwiseHashes(); // ensure seeds set
  auto &ctx = HashContext::instance();

  double h1 = murmur_hash(0, ctx.seed1);
  double h2 = murmur_hash(123456, ctx.seed1);
  double h3 = murmur_hash(-9999, ctx.seed2);

  REQUIRE(h1 >= 0.0);
  REQUIRE(h1 <= 1.0);
  REQUIRE(h2 >= 0.0);
  REQUIRE(h2 <= 1.0);
  REQUIRE(h3 >= 0.0);
  REQUIRE(h3 <= 1.0);
}
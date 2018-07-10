#define TEST_CURVE25519_FIELD_OPS 1

#include <Crypto.h>
#include <BLAKE2s.h>
#include <string.h>
#include <avr/pgmspace.h>

//This library is just used for fast arithmetic operations, calculations are all on FOURQ CURVE.
#include <Curve25519.h>
#include <BigNumberUtil.h>
#include <utility/LimbUtil.h>


BLAKE2s blake;

//static limb_t const numQ[NUM_LIMBS_256BIT] PROGMEM = {
//    LIMB_PAIR(0x2FB2540E, 0xC7768CE7), LIMB_PAIR(0xDFBD004D, 0xFE0F7999),
//    LIMB_PAIR(0xF0539782, 0x9CBC14E5), LIMB_PAIR(0x0029CBC1, 0x4E5E0A72)
//};

static limb_t const numQ[NUM_LIMBS_256BIT] PROGMEM = {
    LIMB_PAIR(0xC7768CE7, 0x2FB2540E), LIMB_PAIR(0xFE0F7999, 0xDFBD004D),
    LIMB_PAIR(0x9CBC14E5, 0xF0539782), LIMB_PAIR(0x4E5E0A72, 0x0029CBC1)
};

void reduceQ(limb_t *result, limb_t *r)
{
    // Algorithm from: http://en.wikipedia.org/wiki/Barrett_reduction
    //
    // We assume that r is less than or equal to (q - 1)^2.
    //
    // We want to compute result = r mod q.  Find the smallest k such
    // that 2^k > q.  In our case, k = 253.  Then set m = floor(4^k / q)
    // and let r = r - q * floor(m * r / 4^k).  This will be the result
    // or it will be at most one subtraction of q away from the result.
    //
    // Note: 4^k = 4^253 = 2^506 = 2^512/2^6.  We can more easily compute
    // the result we want if we set m = floor(4^k * 2^6 / q) instead and
    // then r = r - q * floor(m * r / 2^512).  Because the slight extra
    // precision in m, r is at most two subtractions of q away from the
    // final result.
//    static limb_t const numM[NUM_LIMBS_256BIT + 1] PROGMEM = {
//        LIMB_PAIR(0x0A2C131B, 0xED9CE5A3), LIMB_PAIR(0x086329A7, 0x2106215D),
//        LIMB_PAIR(0xFFFFFFEB, 0xFFFFFFFF), LIMB_PAIR(0xFFFFFFFF, 0xFFFFFFFF),
//        0x0F
//    };

    static limb_t const numM[NUM_LIMBS_256BIT + 2] PROGMEM = {
        LIMB_PAIR(0x81F6A449, 0xE6858D04), LIMB_PAIR(0x72291EA1, 0x809210C3),
        LIMB_PAIR(0x00002251, 0x00000000), LIMB_PAIR(0x00000000, 0x00000000),
        0x0620
    };
    limb_t temp[NUM_LIMBS_512BIT + NUM_LIMBS_256BIT + 1];

    // Multiply r by m.
    BigNumberUtil::mul_P(temp, r, NUM_LIMBS_512BIT, numM, NUM_LIMBS_256BIT + 1);

    // Multiply (m * r) / 2^512 by q and subtract it from r.
    // We can ignore the high words of the subtraction result
    // because they will all turn into zero after the subtraction.
    BigNumberUtil::mul_P(temp, temp + NUM_LIMBS_512BIT, NUM_LIMBS_256BIT + 1,
                         numQ, NUM_LIMBS_256BIT);
    BigNumberUtil::sub(r, r, temp, NUM_LIMBS_256BIT);

    // Perform two subtractions of q from the result to reduce it.
    BigNumberUtil::reduceQuick_P(result, r, numQ, NUM_LIMBS_256BIT);
    BigNumberUtil::reduceQuick_P(result, result, numQ, NUM_LIMBS_256BIT);

    // Clean up and exit.
    clean(temp);
}

//void reduceQ(limb_t *result, limb_t *r){
//01880000000000000000000000000000089460248430DC8A47A879A16341207D
//  static limb_t const numM[NUM_LIMBS_256BIT + 1] PROGMEM = {
//      LIMB_PAIR(0x0A2C131B, 0xED9CE5A3), LIMB_PAIR(0x086329A7, 0x2106215D),
//      LIMB_PAIR(0xFFFFFFEB, 0xFFFFFFFF), LIMB_PAIR(0xFFFFFFFF, 0xFFFFFFFF),
//      0x0F
//  };
//
//  
//}

void setup() {
  Serial.begin(115200);
  Serial.print("Testing MRETA Signing\n");
  
  uint8_t y[32] = {0x54, 0xa2, 0xf8, 0x03, 0x1d, 0x18, 0xac, 0x77, 0xd2, 0x53, 0x92, 0xf2, 0x80, 0xb4, 0xb1, 0x2f, 0xac, 0xf1, 0x29, 0x3f, 0x3a, 0xe6, 0x77, 0x7d, 0x74, 0x15, 0x67, 0x91, 0x99, 0x53, 0x69, 0xc5};
  uint8_t msg[32] = {"a"};
  uint8_t msgCount[2] = {0};
  uint8_t hashed[64], zj[32], cj[32], sj[32];
  unsigned long a1, b1, clockcycle;
  limb_t limby[NUM_LIMBS_256BIT], limbrj[NUM_LIMBS_256BIT], limbsj[NUM_LIMBS_256BIT], limbej[NUM_LIMBS_256BIT], limbtemp[NUM_LIMBS_512BIT];
  BigNumberUtil::unpackLE(limby, NUM_LIMBS_256BIT, y, 32);


//  msgCount[0] = 1;
  a1 = micros();
  msgCount[1] = 0;
  blake.reset(y, sizeof(y), 32);
  blake.update(msgCount, 2);
  blake.finalize(limbrj, 32);

  msgCount[1] = 1;
  blake.reset(y, sizeof(y), 32);
  blake.update(msgCount, 2);
  blake.finalize(zj, 32);

  for(int j = 0; j<32; j++){
      cj[j] = zj[j]^msg[j];
  } 

  blake.reset();
  blake.update(cj, 32);
  blake.finalize(limbej, 32);

  Curve25519::mulNoReduce(limbtemp, limbej, limby);
//  BigNumberUtil::mul(limbtemp, limbej, NUM_LIMBS_256BIT, limby, NUM_LIMBS_256BIT);
  reduceQ(limbtemp, limbtemp);

  if(BigNumberUtil::sub(limbsj, limbrj, limbtemp, NUM_LIMBS_256BIT) == 1){
    Serial.print("Borrow!!!");
    BigNumberUtil::add_P(limbsj, limbsj, numQ, NUM_LIMBS_256BIT);
  }
    

  b1 = micros();
  clockcycle = microsecondsToClockCycles(b1-a1);
  Serial.print("MRETA Signing in: "); Serial.println(clockcycle);
  BigNumberUtil::packLE(sj, 32, limbsj, NUM_LIMBS_256BIT);
  Serial.print("sj = {");
  for (unsigned i = 0; i < 32; ++i) { 
     Serial.print("0x"); Serial.print(sj[i], HEX); Serial.print(", ");
  }
  Serial.println("};");

  Serial.print("cj = {");
  for (unsigned i = 0; i < 32; ++i) { 
     Serial.print("0x"); Serial.print(cj[i], HEX); Serial.print(", ");
  }
  Serial.println("};");
  
}

void loop() {
  // put your main code here, to run repeatedly:

}

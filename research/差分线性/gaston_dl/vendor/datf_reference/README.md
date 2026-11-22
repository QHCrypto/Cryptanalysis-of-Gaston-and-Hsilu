This is the source code repository for the paper "Differential-Linear Cryptanalysis from an Algebraic Perspective", including the code for both theoretical (.sage) and experimental (.cpp) biases.

This repository contains the following files:

   - [4rascon.sage](./4rascon.sage), [5rascon.sage](./5rascon.sage)

     the theoretical bias of the DL distinguisher of Ascon (rounds 0-3.5/4.5)

  - [Verify5rascon.cpp](./Verify5rascon.cpp)
    
    the full DL attack on 5-round Ascon (rounds 0-4.5)
  
  - [3rserpent.sage](./3rserpent.sage), [5rserpent.sage](./5rserpent.sage)
    
    the theoretical bias of the DL distinguisher of 3-round Serpent (rounds 3-5) and the conditional DL distinguisher of 5-round Serpent (rounds 1-5)

  - [VerifySerpentDistinguisher.cpp](./VerifySerpentDistinguisher.cpp)

    the experimental bias of the DL distinguisher of 4-round Serpent (rounds 2-5)

  - [199rGrain128a.sage](./199rGrain128a.sage), [Verify199rGrain128a.cpp](./Verify199rGrain128a.cpp)
    
    the theoretical/experimental bias of the CD distinguisher of 199-round Grain128a

  - [224rGrain128AEADv2.sage](./224rGrain128AEADv2.sage), [Verify224rGrain128AEADv2.cpp](./Verify224rGrain128AEADv2.cpp)
    
    the theoretical/experimental bias of the CD distinguisher of 224-round Grain128AEADv2

  - [Delta288.sage](./Delta288.sage), [VerifyDelta288.cpp](./VerifyDelta288.cpp)
    
    the theoretical/experimental bias of the CD distinguisher of 160-round Grain128a ($\Delta b_{288}$)

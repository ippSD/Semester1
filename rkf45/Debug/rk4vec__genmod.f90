        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 06 11:30:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RK4VEC__genmod
          INTERFACE 
            SUBROUTINE RK4VEC(T0,M,U0,DT,F,U)
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: T0
              REAL(KIND=8) :: U0(M)
              REAL(KIND=8) :: DT
              EXTERNAL F
              REAL(KIND=8) :: U(M)
            END SUBROUTINE RK4VEC
          END INTERFACE 
        END MODULE RK4VEC__genmod

        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 06 11:30:18 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RK4__genmod
          INTERFACE 
            SUBROUTINE RK4(T0,U0,DT,F,U)
              REAL(KIND=8) :: T0
              REAL(KIND=8) :: U0
              REAL(KIND=8) :: DT
              EXTERNAL F
              REAL(KIND=8) :: U
            END SUBROUTINE RK4
          END INTERFACE 
        END MODULE RK4__genmod

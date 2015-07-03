!-----------------------------------------------------------------------
!
!  Hawkmatix Statistics Module for Fortran
!  Official project page: http://www.hawkmatix.com/statistics.html
!
!  Copyright 2014, 2015 Andrew C. Hawkins (andrew.hawkins@hawkmatix.com)
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as
!  published by the Free Software Foundation, either version 3 of the
!  License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this program.  If not, see <http://www.gnu.org/
!  licenses/>.
!
!-----------------------------------------------------------------------

MODULE STAT

IMPLICIT NONE
CONTAINS

!-----------------------------------------------------------------------
!
!  Module Functions
!
!-----------------------------------------------------------------------

FUNCTION ARITHMETIC_MEAN(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the arithmetic mean of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      arithmetic mean will be found.
    !
    !  Retruns:
    !    DOUBLE PRECISION: The arithmetic mean.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: ARITHMETIC_MEAN

    ARITHMETIC_MEAN = SUM(DATASET) / SIZE(DATASET)
END FUNCTION ARITHMETIC_MEAN

FUNCTION GEOMETRIC_MEAN(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the geometric mean of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      geometric mean will be found.
    !
    !  Retruns:
    !    DOUBLE PRECISION: The geometric mean.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: GEOMETRIC_MEAN

    GEOMETRIC_MEAN = PRODUCT(DATASET) ** (1.0D0 / SIZE(DATASET))
END FUNCTION GEOMETRIC_MEAN

FUNCTION HARMONIC_MEAN(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the harmonic mean of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      harmonic mean will be found.
    !
    !  Retruns:
    !    DOUBLE PRECISION: The harmonic mean.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: HARMONIC_MEAN

    HARMONIC_MEAN = SIZE(DATASET) / SUM(1.0D0 / DATASET)
END FUNCTION HARMONIC_MEAN

FUNCTION QUADRATIC_MEAN(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the quadratic mean of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      quadtratic mean will be found.
    !
    !  Retruns:
    !    DOUBLE PRECISION: The quadratic mean.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: QUADRATIC_MEAN

    QUADRATIC_MEAN = DSQRT(SUM(DATASET * DATASET) / SIZE(DATASET))
END FUNCTION QUADRATIC_MEAN

FUNCTION GENERALIZED_MEAN(DATASET, POWER)
    !-------------------------------------------------------------------
    !
    !  Calculate the generalized mean of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      generalized mean will be found.
    !
    !  Retruns:
    !    DOUBLE PRECISION: The generalized mean.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: POWER
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: GENERALIZED_MEAN

    GENERALIZED_MEAN = (SUM(DATASET ** POWER) / SIZE(DATASET)) ** &
        (1.0D0 / POWER)
END FUNCTION GENERALIZED_MEAN

FUNCTION WEIGHTED_MEAN(DATASET, WEIGHT)
    !-------------------------------------------------------------------
    !
    !  Calculate the weighted mean of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      weighted mean will be found.
    !    WEIGHT (DOUBLE PRECISION):
    !
    !  Retruns:
    !    DOUBLE PRECISION: The weighted mean.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET, WEIGHT
    DOUBLE PRECISION :: WEIGHTED_MEAN

    WEIGHTED_MEAN = SUM(DATASET * WEIGHT) / SUM(WEIGHT)
END FUNCTION WEIGHTED_MEAN

FUNCTION MEDIAN(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the median of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the median
    !      will be found.
    !
    !  Retruns:
    !    DOUBLE PRECISION: The median.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: MEDIAN
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SORTED
    INTEGER :: N, MP

    N = SIZE(DATASET)
    ALLOCATE(SORTED(N))
    SORTED = SORT(N, DATASET)

    IF (MOD(N, 2) == 0) THEN
        MP = N / 2
        MEDIAN = (SORTED(MP) + SORTED(MP + 1)) / 2.0D0
    ELSE
        MP = N / 2
        MEDIAN = SORTED(MP + 1)
    END IF
END FUNCTION MEDIAN

FUNCTION MODE(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the mode of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the mode
    !      will be found.
    !
    !  Retruns:
    !    DOUBLE PRECISION: The mode.
    !
    !  To Do:
    !    Find the modes of multimodal datasets.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: MODE
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SORTED, COUNT
    INTEGER :: N, i

    N = SIZE(DATASET)
    ALLOCATE(SORTED(N), COUNT(N))
    SORTED = SORT(N, DATASET)

    COUNT(1) = 1

    DO i = 1, N - 1
        IF (SORTED(i) == SORTED(i + 1)) THEN
            COUNT(i + 1) = COUNT(i) + 1
        ELSE
            COUNT(i + 1) = 1
        END IF
    END DO

    MODE = SORTED(MAX_LOC(N, COUNT))
END FUNCTION MODE

FUNCTION QUARTILE(DATASET, Q)
    !-------------------------------------------------------------------
    !
    !  Calculate the mode of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      quartile will be found.
    !    Q (INTEGER): The quartile to find (between 1 and 3).
    !
    !  Retruns:
    !    DOUBLE PRECISION: The quartile.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    INTEGER :: Q
    DOUBLE PRECISION :: QUARTILE
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SORTED
    INTEGER :: N, i

    N = Size(DATASET)
    i = Q * (N - 1)
    ALLOCATE(SORTED(N))
    SORTED = SORT(N, DATASET)

    QUARTILE = SORTED(i / 4 + 1) + (MOD(i, 4) / 4) * &
        (SORTED(i / 4 + 2) - SORTED(i / 4 + 1))
END FUNCTION QUARTILE

FUNCTION SRANGE(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the span range of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the span
    !      range will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The span range.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: SRANGE

    SRANGE = MAXVAL(DATASET) - MINVAL(DATASET)
END FUNCTION SRANGE

FUNCTION MIDRANGE(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the midrange of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      midrange will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The midrange.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: MIDRANGE

    MIDRANGE = (MAXVAL(DATASET) + MINVAL(DATASET)) / 2.0D0
END FUNCTION MIDRANGE

FUNCTION IQRANGE(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the interquartile range of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      interquartile range will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The interquartile range.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: IQRANGE
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SORTED
    INTEGER :: N
    DOUBLE PRECISION :: LOWER, UPPER

    N = SIZE(DATASET)
    ALLOCATE(SORTED(N))
    SORTED = SORT(N, DATASET)

    IF (MOD(N, 2) == 0) THEN
        LOWER = MEDIAN(SORTED(1 : N / 2))
        UPPER = MEDIAN(SORTED(N / 2 + 1 : N))
    ELSE
        LOWER = MEDIAN(SORTED(1 : N / 2 + 1))
        UPPER = MEDIAN(SORTED(N / 2 + 1 : N))
    END IF

    IQRANGE = UPPER - LOWER
END FUNCTION IQRANGE

FUNCTION TRIMEAN(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the trimean of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      trimean will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The trimean.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: TRIMEAN

    TRIMEAN = (QUARTILE(DATASET, 1) + 2 * QUARTILE(DATASET, 2) + &
        QUARTILE(DATASET, 3)) / 4.0D0
END FUNCTION TRIMEAN

FUNCTION VARIANCE(DATASET, POPULATION)
    !-------------------------------------------------------------------
    !
    !  Calculate the variance of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      variance will be found.
    !    POPULATION (LOGICAL): True if the dataset represents a
    !      population, false otherwise.
    !
    !  Returns:
    !    DOUBLE PRECISION: The variance.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: POPULATION
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: VARIANCE
    INTEGER :: i, length
    DOUBLE PRECISION :: bulk, mean

    length = SIZE(DATASET)
    bulk = 0.0D0
    mean = ARITHMETIC_MEAN(DATASET)

    DO i = 1, length
        bulk = bulk + (DATASET(i) - mean) ** 2
    END DO

    IF (POPULATION .EQV. .FALSE.) THEN
        VARIANCE = bulk / (length - 1)
    ELSE
        VARIANCE = bulk / length
    END IF
END FUNCTION VARIANCE

FUNCTION STANDARD_DEVIATION(DATASET, POPULATION)
    !-------------------------------------------------------------------
    !
    !  Calculate the standard deviation of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      standard deviation will be found.
    !    POPULATION (LOGICAL): True if the dataset represents a
    !      population, false otherwise.
    !
    !  Returns:
    !    DOUBLE PRECISION: The standard deviation.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: POPULATION
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: STANDARD_DEVIATION

    STANDARD_DEVIATION = DSQRT(VARIANCE(DATASET, POPULATION))
END FUNCTION STANDARD_DEVIATION

FUNCTION MEAN_ABSOLUTE_DEVIATION(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the mean absolute deviation of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the mean
    !      absolute deviation will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The mean absolute deviation.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: MEAN_ABSOLUTE_DEVIATION
    INTEGER :: i, length
    DOUBLE PRECISION :: bulk, mean

    length = SIZE(DATASET)
    bulk = 0.0D0
    mean = ARITHMETIC_MEAN(DATASET)

    DO i = 1, length
        bulk = bulk + DABS(DATASET(i) - mean)
    END DO

    MEAN_ABSOLUTE_DEVIATION = bulk / length
END FUNCTION MEAN_ABSOLUTE_DEVIATION

FUNCTION MEAN_DIFFERENCE(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the mean difference of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the span
    !      range will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The span range.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: MEAN_DIFFERENCE
    INTEGER :: i, j, length
    DOUBLE PRECISION :: bulk

    length = SIZE(DATASET)
    bulk = 0.0D0

    DO i = 1, length
        DO j = 1, length
            bulk = bulk + DABS(DATASET(i) - DATASET(j))
        END DO
    END DO

    MEAN_DIFFERENCE = bulk / (length ** 2)
END FUNCTION MEAN_DIFFERENCE

FUNCTION RELATIVE_MEAN_DIFFERENCE(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the relative mean difference of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      relative mean difference will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The relative mean difference.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: RELATIVE_MEAN_DIFFERENCE

    RELATIVE_MEAN_DIFFERENCE = MEAN_DIFFERENCE(DATASET) / &
        ARITHMETIC_MEAN(DATASET)
END FUNCTION RELATIVE_MEAN_DIFFERENCE

FUNCTION MEDIAN_ABSOLUTE_DEVIATION(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the median absolute deviation of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      median absolute deviation will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The median absolute deviation.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: MEDIAN_ABSOLUTE_DEVIATION, MED
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DIFF
    INTEGER :: N, i

    N = SIZE(DATASET)
    ALLOCATE(DIFF(N))
    MED = MEDIAN(DATASET)

    DO i = 1, N
        DIFF(i) = DABS(DATASET(i) - MED)
    END DO

    MEDIAN_ABSOLUTE_DEVIATION = MEDIAN(DIFF)
END FUNCTION MEDIAN_ABSOLUTE_DEVIATION

FUNCTION SKEWNESS(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the skewness of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      skewness will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The skewness.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: SKEWNESS
    INTEGER :: i, length
    DOUBLE PRECISION :: bulk, mean

    length = SIZE(DATASET)
    bulk = 0.0D0
    mean = ARITHMETIC_MEAN(DATASET)

    DO i = 1, length
        bulk = bulk + (DATASET(i) - mean) ** 3
    END DO

    SKEWNESS = (bulk / length) / VARIANCE(DATASET, .FALSE.) ** &
        (3 / 2.0D0)
END FUNCTION SKEWNESS

FUNCTION KURTOSIS(DATASET)
    !-------------------------------------------------------------------
    !
    !  Calculate the kurtosis of a dataset.
    !
    !  Args:
    !    DATASET (DOUBLE PRECISION(:)): The dataset for which the
    !      kurtosis will be found.
    !
    !  Returns:
    !    DOUBLE PRECISION: The kurtosis.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DATASET
    DOUBLE PRECISION :: KURTOSIS
    INTEGER :: i, length
    DOUBLE PRECISION :: bulk, mean

    length = SIZE(DATASET)
    bulk = 0.0D0
    mean = ARITHMETIC_MEAN(DATASET)

    DO i = 1, length
        bulk = bulk + (DATASET(i) - mean) ** 4
    END DO

    KURTOSIS = (bulk / length) / VARIANCE(DATASET, .FALSE.) ** 2
END FUNCTION KURTOSIS

FUNCTION FACTORIAL(N)
    !-------------------------------------------------------------------
    !
    !  Calculate the factorial of a number.
    !
    !  Args:
    !    N (INTEGER): The factorial number.
    !
    !  Returns:
    !    DOUBLE PRECISION: The factorial.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION :: FACTORIAL

    FACTORIAL = GAMMA(N + 1.0D0)
END FUNCTION FACTORIAL

FUNCTION PFACTORIAL(N, M)
    !-------------------------------------------------------------------
    !
    !  Calculate the partial factorial of a number.
    !
    !  Args:
    !    N (INTEGER): The factorial number.
    !    M (INTEGER): The lower bound of the factorial.
    !
    !  Returns:
    !    DOUBLE PRECISION: The factorial.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, M
    DOUBLE PRECISION :: PFACTORIAL
    INTEGER :: i
    DOUBLE PRECISION :: fact

    fact = 1.0D0

    DO i = N, M, -1
        fact = fact * i
    END DO

    PFACTORIAL = fact
END FUNCTION PFACTORIAL

FUNCTION PERMUTATION(N, K)
    !-------------------------------------------------------------------
    !
    !  Calculate the permutation.
    !
    !  Args:
    !    N (INTEGER): The N parameter.
    !    K (INTEGER): The K parameter.
    !
    !  Returns:
    !    DOUBLE PRECISION: The permutation.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, K
    DOUBLE PRECISION :: PERMUTATION

    IF (N - K > 10) THEN
        PERMUTATION = FACTORIAL(N) / FACTORIAL(N - K)
    ELSE
        PERMUTATION = PFACTORIAL(N, N - K + 1)
    END IF
END FUNCTION PERMUTATION

FUNCTION COMBINATION(N, K)
    !-------------------------------------------------------------------
    !
    !  Calculate the combination.
    !
    !  Args:
    !    N (INTEGER): The N parameter.
    !    K (INTEGER): The K parameter.
    !
    !  Returns:
    !    DOUBLE PRECISION: The combination.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, K
    DOUBLE PRECISION :: COMBINATION

    IF (N - K > 10) THEN
        COMBINATION = FACTORIAL(N) / (FACTORIAL(K) * FACTORIAL(N - K))
    ELSE
        IF (K > N - K) THEN
            COMBINATION = PFACTORIAL(N, K + 1) / FACTORIAL(N - K)
        ELSE
            COMBINATION = PFACTORIAL(N, N - K + 1) / FACTORIAL(K)
        END IF
    END IF
END FUNCTION COMBINATION

!-----------------------------------------------------------------------
!
!  Helper Functions
!
!-----------------------------------------------------------------------

FUNCTION SORT(N, DATASET) RESULT(SORTED)
    !-------------------------------------------------------------------
    !
    !  Sort a one-dimensional array.
    !
    !  Args:
    !    N (INTEGER): The size of the dataset.
    !    DATASET (DOUBLE PRECISION(N)): The dataset to be sorted.
    !
    !  Returns:
    !    DOUBLE PRECISION(N): The sorted dataset.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: DATASET
    DOUBLE PRECISION, DIMENSION(N) :: SORTED
    INTEGER :: i, j, a

    SORTED = DATASET

    DO j = 2, N
        a = SORTED(j)
        DO i = j - 1, 1, -1
            IF (SORTED(i) <= a) GOTO 10
            SORTED(i + 1) = SORTED(i)
        END DO
        i = 0
10      SORTED(i + 1) = a
    END DO
END FUNCTION SORT

FUNCTION MAX_LOC(N, DATASET)
    !-------------------------------------------------------------------
    !
    !  Determine the maximum value and location of a dataset.
    !
    !  Args:
    !    N (INTEGER): The size of the dataset.
    !    DATASET (DOUBLE PRECISION(N)): The dataset to find the maximum
    !      of.
    !
    !  Returns:
    !    DOUBLE PRECISION: The first maximum value's index.
    !
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: DATASET
    INTEGER :: MAX_LOC
    DOUBLE PRECISION :: a
    INTEGER :: b, i

    a = -HUGE(0.0D0)
    MAX_LOC = 0
    DO i = 1, N
        IF (DATASET(i) > a) THEN
            a = DATASET(i)
            MAX_LOC = i
        END IF
    END DO
END FUNCTION MAX_LOC

END MODULE STAT

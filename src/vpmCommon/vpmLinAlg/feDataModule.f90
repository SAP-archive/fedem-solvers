!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module FEData

  !!============================================================================
  !! This module defines a data type that is used to transport the FE data
  !! through Fedem into the DNVS _Linear Algebra Modules_ subroutines.
  !! Note: The module *must* be named FEData and the data type FEDataInput,
  !! as these names are also used within the DNVS library. However, the contents
  !! of the data type is arbitrary and can be application (Fedem) specific.
  !!============================================================================

  use SamModule, only : SamType

  implicit none

  type FEDataInput
     type(SamType), pointer :: sam => null()
  end type FEDataInput

end module FEData

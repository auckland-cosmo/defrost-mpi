! Define the model here...

		
! header for models. Change based on the model. Make sure new parameters are specified in the parameter file.

subroutine modelsummary(fd)
        integer fd

        write (fd,'(a,3(f0.5,a))') "# V(phi,psi) = ", &
                sqrt(m2phi), "^2*(phi^2) this is garbage"        
!                sqrt(m2psi), "^2*psi^2/2 + ",         &
!                sqrt(g2), "^2*phi^2*psi^2/2 + ",	  &
!				lambda, "*phi^4/4"
				
end subroutine modelsummary

! calculates the value of the potential at each point

function modelv(hr, i, j, k)
        real modelv
        real, dimension(fields,sistart(1):siend(1),sistart(2):siend(2),sistart(3):siend(3)) :: hr; integer i, j, k
		 modelv = (hr(phi,i,j,k)**2-(lambdapf*ll/2.0)*(mpl**2/m2phi)*hr(phi,i,j,k)**4+(gpf*gg**2/3.0)*(mpl**4/m2phi**2)*hr(phi,i,j,k)**6) ! This is 2*V. a factor of mpl**2 is ignored
!        modelv = (m2phi + g2*hr(psi,i,j,k)**2)*hr(phi,i,j,k)**2 + m2psi*hr(psi,i,j,k)**2 + lambda*hr(phi,i,j,k)**4/4  ! change this line for each model
end function modelv

! calculates the potential derivatives at each point wrt the fields. modeldv is a 2 by 1 array. First row is (dV/d\phi)/phi and second row is (dV/dpsi)/psi

function modeldv(hr, i, j, k)
        real, dimension(fields) :: modeldv
        real, dimension(fields,sistart(1):siend(1),sistart(2):siend(2),sistart(3):siend(3)) :: hr; integer i, j, k

        ! basis vectors for vectorizing potential derivative
        real, dimension(fields), parameter :: V1 = 1
!(/1,0/), V2 = (/0,1/)
			
	 modeldv = (1 -  lambdapf*ll*(mpl**2/m2phi)*hr(phi,i,j,k)**2 +   (gpf*gg**2)*(mpl**4/m2phi**2)*hr(phi,i,j,k)**4)*V1 
!        modeldv = (1 + g2*hr(psi,i,j,k)**2 + lambda*hr(phi,i,j,k)**2)*V1 + (m2psi + g2*hr(phi,i,j,k)**2)*V2 
end function modeldv

! Contribution to the effective mass from the model
function modelm2eff()
        real, dimension(fields) :: modelm2eff
        modelm2eff = 0.
!(/ 0.0, g2*phi0**2 /)
	return
end function modelm2eff


! hetero with no intercepts

subroutine pmqmlam(x,betax,residualsx,n,p,k,tau,m,a,eps,maxin,maxout,lambdaset,pentype,ll)
implicit none
integer*4 :: n, p, K, m, maxin, maxout, i, pentype, ll
real*8 :: x(n,p*K), betax(ll, p*K*m), residualsx(n*ll,m*K), a, eps, &
  tau(m), lambdaset(ll), beta(p*K*m), residuals(n, m*K), lambda

  do i = 1, ll
  
    lambda = lambdaset(i)
    beta = betax(i, :)
    residuals = residualsx( (n * (i - 1) + 1) : (n * i), : )
    
    call multicqpenheteronointer(x,beta,residuals,n,p,k,tau,m,a,eps,maxin,maxout,lambda,pentype)
    
    betax(i, :) = beta
  
  end do

return
end subroutine pmqmlam


subroutine multicqpenheteronointer(x,beta,residuals,n,p,k,tau,m,a,eps,maxin,maxout,lambda,pentype)
implicit none
integer*4 :: n, p, K, m, maxin, count, iter, rep, maxout, i, pentype, rep1, count1
real*8 :: x(n,p*K), beta(p*K*m), beta0(p*K*m), penweight(p), residuals(n,m*K), a, eps, &
  xtemp(n,p),betatemp(p), residualstemp(n), distance, tau(m), lambda, tau1


i=1
distance=eps+1.d0

do while((i<=maxout) .AND. (distance>=eps))
  
  beta0=beta

  do iter=1,p

    penweight(iter)=0.d0

    do rep=1,K
    
      do rep1=1,m

        penweight(iter)=penweight(iter)+abs(beta((rep-1)*m*p+(rep1-1)*p+iter))
        ! penweight(iter)=penweight(iter)+beta((rep-1)*p+iter)**2
      
      end do
      
    end do

  end do

  ! penweight=sqrt(penweight)

  call pender(penweight,p,a,lambda,pentype)
  penweight=penweight*real(n,8)

  do count=1,K
    
    do count1=1,m

      xtemp=x(:,(p*(count-1)+1):(p*count))
      betatemp=beta((p*m*(count-1)+p*(count1-1)+1):(p*m*(count-1)+p*count1))
      residualstemp=residuals(:,(m*(count-1)+count1))
      tau1=tau(count1)

      call qpennointer(xtemp,betatemp,penweight,residualstemp,n,p,tau1,eps,maxin)

      beta((p*m*(count-1)+p*(count1-1)+1):(p*m*(count-1)+p*count1))=betatemp
      residuals(:,(m*(count-1)+count1))=residualstemp
      
    end do

  end do

  distance=sqrt(sum((beta-beta0)**2))

  i=i+1

end do

if((i==maxout) .AND. (distance>eps)) then

  print*,"bad"

end if

return
end subroutine multicqpenheteronointer



subroutine qpennointer(x,beta,penweight,residuals,n,p,tau,eps,maxin)
implicit none
integer*4 :: n, p, m, maxin, iter, col, row, nonzerox, count3, count4
real*8 :: x(n,p), beta(p), penweight(p), residuals(n), tau, eps, betavecDiff, &
  weight_sum, weight_sum1, weight_sum2, pre_value_final(n+1), weight(n+1), temp1, temp2, newBeta, &
  betaDiff


do iter=1,maxin
  col=0
  betavecDiff=0.d0

  do col=1,p

    weight_sum=0.d0
    nonzerox=0

    do row=1,n

      if (x(row,col)/=0.d0) then
        nonzerox=nonzerox+1
        pre_value_final(nonzerox)=residuals(row)

        if (pre_value_final(nonzerox)>0.d0) then
          weight(nonzerox)=abs(x(row,col))*(tau)
        else
          weight(nonzerox)=abs(x(row,col))*(1-tau)
        end if

        pre_value_final(nonzerox)=(pre_value_final(nonzerox)+x(row,col)*beta(col))/(x(row,col))

        weight_sum=weight_sum+weight(nonzerox)
      end if

    end do

    nonzerox=nonzerox+1
    pre_value_final(nonzerox)=0.d0

    weight(nonzerox)=penweight(col)
    weight_sum=weight_sum+weight(nonzerox)

    call xssort(pre_value_final,weight,nonzerox,2)

    count3=0
    temp1=0.d0
    weight_sum=weight_sum/2.d0
    do while(temp1<=weight_sum )
  
      count3=count3+1
      temp1=temp1+weight(count3)

    end do	

    count4=nonzerox+1
    temp2=0.d0
    do while(temp2<=weight_sum )
  
      count4=count4-1
    
      temp2=temp2+weight(count4)

    end do

    newBeta=(pre_value_final(count3)+pre_value_final(count4))/2.d0

    betaDiff=beta(col)-newBeta
    betavecDiff=betavecDiff+betaDiff*betaDiff
    beta(col)=newBeta

    if (betaDiff/=0.d0) then

      do row=1,n

        residuals(row)=residuals(row)+x(row,col)*betaDiff

      end do

    end if

  end do

  if(betavecDiff<eps*eps) then

    exit

  end if

end do

return
end subroutine qpennointer


! hetero with intercepts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cross validation. 
subroutine cvpqr(x, y,loss,betaset,intvalset,residuals,n,p,tau,a,eps,maxin,maxout, &
                lambdaset,pentype,ll,nfold,n1,n2)
implicit none
integer*4 :: n, p, K, m, maxin, maxout, i, pentype, ll, nfold, n1, n2, &
             ind1(n-n1), z, j, r
real*8 :: x(n,p), betax(ll, p), a, eps, tau(1), lambdaset(ll), &
          betaset(ll, p), residuals(n, 1), intvalset(ll,1), & 
          intvalx(ll,1), xx1(n-n1,p), xx2(n-n2,p), residuals1(n-n1, 1), &
          residuals2(n-n2, 1), loss(ll), y(n), xv1(n1,p), xv2(n2,p), &
          yv1(n1, 1), yv2(n2, 1), error1(n1, 1), error2(n2, 1), slope(p, 1), indicator1(n1), &
          indicator2(n2)
  
  loss = 0.d0
  m = 1
  K = 1
  
  xx1 = x((n1 + 1) : n, :)
  residuals1 = residuals((n1 + 1) : n, :)
  betax = betaset
  intvalx = intvalset
  
  call pqrlambda(xx1,betax,intvalx,residuals1,n-n1,p,k,tau,m,a,eps,maxin,maxout, &
                      lambdaset,pentype,ll)
  yv1(:, 1) = y(1 : n1)
  xv1 = x(1 : n1, :)
                      
  do j = 1, ll
  
    slope(:, 1) = betax(j, :)
  
    error1 = yv1 - matmul( xv1, slope ) - intvalx(j, 1)
    
    do r = 1, n1
      
      indicator1(r) = merge(1.d0, 0.d0, error1(r, 1) < 0.d0)
    
    end do
    
    loss(j) = loss(j) + sum(error1(:, 1) * (tau(1) - indicator1))
  
  end do
                      
  
  do i = 2, (nfold - 1)
  
    ind1 = (/ (z, z = 1, (i - 1) * n1), (z, z = i * n1 + 1, n)/)
    
    xx1 = x(ind1, :)
    residuals1 = residuals(ind1, :)
    betax = betaset
    intvalx = intvalset
  
    call pqrlambda(xx1,betax,intvalx,residuals1,n-n1,p,k,tau,m,a,eps,maxin,maxout, &
                      lambdaset,pentype,ll)
                      
    yv1(:, 1) = y( ((i - 1) * n1 + 1) : (i * n1) )
    xv1 = x( ((i - 1) * n1 + 1) : (i * n1), : )
                      
    do j = 1, ll
  
      slope(:, 1) = betax(j, :)
  
      error1 = yv1 - matmul( xv1, slope ) - intvalx(j, 1)
  
      do r = 1, n1
      
        indicator1(r) = merge(1.d0, 0.d0, error1(r, 1) < 0.d0)
    
      end do
    
      loss(j) = loss(j) + sum(error1(:, 1) * (tau(1) - indicator1))
  
    end do  
    
  end do
  
  xx2 = x(1 : (n - n2), :)
  residuals2 = residuals(1 : (n - n2), :)
  betax = betaset
  intvalx = intvalset
  
  call pqrlambda(xx2,betax,intvalx,residuals2,n-n2,p,k,tau,m,a,eps,maxin,maxout, &
                      lambdaset,pentype,ll)
                      
  yv2(:, 1) = y( (n - n2 + 1) : n )
  xv2 = x( (n - n2 + 1) : n, : )
                      
  do j = 1, ll
  
    slope(:, 1) = betax(j, :)
  
    error2 = yv2 - matmul( xv2, slope ) - intvalx(j, 1)
  
    do r = 1, n2
      
      indicator2(r) = merge(1.d0, 0.d0, error2(r, 1) < 0.d0)
    
    end do
    
    loss(j) = loss(j) + sum(error2(:, 1) * (tau(1) - indicator2))
  
  end do  

return
end subroutine cvpqr

subroutine pqrlambda(x,betax,intvalx,residualsx,n,p,k,tau,m,a,eps,maxin,maxout, &
                      lambdaset,pentype,ll)
implicit none
integer*4 :: n, p, K, m, maxin, maxout, i, pentype, ll, continue
real*8 :: x(n,p*K), betax(ll, p*K*m), a, eps, residualsx(n, m*K), &
  tau(m), lambdaset(ll), beta(p*K*m), residuals(n, m*K), lambda, intval(m*K), intvalx(ll,m*K)
  
  continue = 1
  
  do i = 1, ll
    
    if(continue == 1) then
    
      lambda = lambdaset(i)
      beta = betax(i, :)
      intval = intvalx(i, :)
      residuals = residualsx
    
      call multicqpenhetero(x,beta,intval,residuals,n,p,k,tau,m,a,eps,maxin,maxout,lambda,pentype)
    
      betax(i, :) = beta
      intvalx(i, :) = intval
      
    else
    
      betax(i, :) = betax(i - 1, :)
      intvalx(i, :) = intvalx(i - 1, :)
    
    end if
    
    if(all(beta == 0)) then

      continue = 0
    
    else
    
      continue = 1

    end if
    
  
  end do

return
end subroutine pqrlambda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine pmqmlamint(x,betax,intvalx,residualsx,n,p,k,tau,m,a,eps,maxin,maxout, &
                      lambdaset,pentype,ll)
implicit none
integer*4 :: n, p, K, m, maxin, maxout, i, pentype, ll, continue
real*8 :: x(n,p*K), betax(ll, p*K*m), residualsx(n*ll,m*K), a, eps, &
  tau(m), lambdaset(ll), beta(p*K*m), residuals(n, m*K), lambda, intval(m*K), intvalx(ll,m*K)
  
  continue = 1
  
  do i = 1, ll
    
    if(continue == 1) then
    
      lambda = lambdaset(i)
      beta = betax(i, :)
      intval = intvalx(i, :)
      residuals = residualsx( (n * (i - 1) + 1) : (n * i), : )
    
      call multicqpenhetero(x,beta,intval,residuals,n,p,k,tau,m,a,eps,maxin,maxout,lambda,pentype)
    
      betax(i, :) = beta
      intvalx(i, :) = intval
      
    else
    
      betax(i, :) = betax(i - 1, :)
      intvalx(i, :) = intvalx(i - 1, :)
    
    end if
    
    if(all(beta == 0)) then

      continue = 0
    
    else
    
      continue = 1

    end if
    
  
  end do

return
end subroutine pmqmlamint


subroutine multicqpenhetero(x,beta,intval,residuals,n,p,k,tau,m,a,eps,maxin,maxout,lambda,pentype)
implicit none
integer*4 :: n, p, K, m, maxin, count, iter, rep, maxout, i, pentype, rep1, count1
real*8 :: x(n,p*K), beta(p*K*m), intval(m*K), beta0(p*K*m), intval0(m*K), penweight(p), residuals(n,m*K), a, eps, &
  xtemp(n,p),betatemp(p),intvaltemp, residualstemp(n), distance, tau(m), lambda, tau1


i=1
distance=eps+1.d0

do while((i<=maxout) .AND. (distance>=eps))
  
  beta0=beta
  intval0=intval

  do iter=1,p

    penweight(iter)=0.d0

    do rep=1,K
    
      do rep1=1,m

        penweight(iter)=penweight(iter)+abs(beta((rep-1)*m*p+(rep1-1)*p+iter))
        ! penweight(iter)=penweight(iter)+beta((rep-1)*p+iter)**2
      
      end do
      
    end do

  end do

  ! penweight=sqrt(penweight)

  call pender(penweight,p,a,lambda,pentype)
  penweight=penweight*real(n,8)

  do count=1,K
    
    do count1=1,m

      xtemp=x(:,(p*(count-1)+1):(p*count))
      betatemp=beta((p*m*(count-1)+p*(count1-1)+1):(p*m*(count-1)+p*count1))
      intvaltemp=intval((count-1)*m+count1)
      residualstemp=residuals(:,(m*(count-1)+count1))
      tau1=tau(count1)

      call qpen(xtemp,betatemp,intvaltemp,penweight,residualstemp,n,p,tau1,eps,maxin)

      beta((p*m*(count-1)+p*(count1-1)+1):(p*m*(count-1)+p*count1))=betatemp
      intval((count-1)*m+count1)=intvaltemp
      residuals(:,(m*(count-1)+count1))=residualstemp
      
    end do

  end do

  distance=sqrt(sum((beta-beta0)**2)+sum((intval0-intval)**2))

  i=i+1

end do

if((i==maxout) .AND. (distance>eps)) then

  print*,"bad"

end if

return
end subroutine multicqpenhetero



subroutine qpen(x,beta,intval,penweight,residuals,n,p,tau,eps,maxin)
implicit none
integer*4 :: n, p, m, maxin, iter, col, row, nonzerox, count3, count4
real*8 :: x(n,p), beta(p), intval, penweight(p), residuals(n), tau, eps, betavecDiff, &
  weight_sum, weight_sum1, weight_sum2, pre_value_final(n+1), weight(n+1), temp1, temp2, newBeta, &
  betaDiff, pre_int_final(n), weight_int(n)


do iter=1,maxin
  col=0
  betavecDiff=0.d0

  do col=1,p

    weight_sum=0.d0
    nonzerox=0

    do row=1,n

      if (x(row,col)/=0.d0) then
        nonzerox=nonzerox+1
        pre_value_final(nonzerox)=residuals(row)

        if (pre_value_final(nonzerox)>0.d0) then
          weight(nonzerox)=abs(x(row,col))*(tau)
        else
          weight(nonzerox)=abs(x(row,col))*(1-tau)
        end if

        pre_value_final(nonzerox)=(pre_value_final(nonzerox)+x(row,col)*beta(col))/(x(row,col))

        weight_sum=weight_sum+weight(nonzerox)
      end if

    end do

    nonzerox=nonzerox+1
    pre_value_final(nonzerox)=0.d0

    weight(nonzerox)=penweight(col)
    weight_sum=weight_sum+weight(nonzerox)

    call xssort(pre_value_final,weight,nonzerox,2)

    count3=0
    temp1=0.d0
    weight_sum=weight_sum/2.d0
    do while(temp1<=weight_sum )
  
      count3=count3+1
      temp1=temp1+weight(count3)

    end do	

    count4=nonzerox+1
    temp2=0.d0
    do while(temp2<=weight_sum )
  
      count4=count4-1
    
      temp2=temp2+weight(count4)

    end do

    newBeta=(pre_value_final(count3)+pre_value_final(count4))/2.d0

    betaDiff=beta(col)-newBeta
    betavecDiff=betavecDiff+betaDiff*betaDiff
    beta(col)=newBeta

    if (betaDiff/=0.d0) then

      do row=1,n

        residuals(row)=residuals(row)+x(row,col)*betaDiff

      end do

    end if

  end do

  ! intercepts

  weight_sum=0.d0

  do row=1,n

    pre_int_final(row)=residuals(row)+intval
    weight_int(row)=1.d0
    weight_sum=weight_sum+weight_int(row)

  end do

  call xssort(pre_int_final,weight_int,n,2)

  count3=0
  temp1=0.d0
  weight_sum1=weight_sum*tau
  do while(temp1<=weight_sum1 )
  
      count3=count3+1
      temp1=temp1+weight_int(count3)

  end do	

  count4=n+1
  temp2=0.d0
  weight_sum2=weight_sum*(1.d0-tau)
  do while(temp2<=weight_sum2 )
  
    count4=count4-1
    temp2=temp2+weight_int(count4)

  end do

  newBeta=(pre_int_final(count3)+pre_int_final(count4))/2.d0

  betaDiff=intval-newBeta
  betavecDiff=betavecDiff+betaDiff*betaDiff
  intval=newBeta

  if (betaDiff/=0.d0) then

    do row=1,n

      residuals(row)=residuals(row)+betaDiff

    end do
 
  end if

  if(betavecDiff<eps*eps) then

    exit

  end if

end do


return
end subroutine qpen











! homo
subroutine multicqpen(x,beta,intval,residuals,n,p,k,tau,m,a,eps,maxin,maxout,lambda,pentype)
implicit none
integer*4 :: n, p, K, m, maxin, count, iter, rep, maxout, i, pentype
real*8 :: x(n,p*K), beta(p*K), intval(m*K), beta0(p*K), intval0(m*K), penweight(p), residuals(n,m*K), a, eps, &
  xtemp(n,p),betatemp(p),intvaltemp(m), residualstemp(n,m), distance, tau(m), lambda


i=1
distance=eps+1.d0

do while((i<=maxout) .AND. (distance>=eps))
  
  beta0=beta
  intval0=intval

  do iter=1,p

    penweight(iter)=0.d0

    do rep=1,K

      penweight(iter)=penweight(iter)+abs(beta((rep-1)*p+iter))
      ! penweight(iter)=penweight(iter)+beta((rep-1)*p+iter)**2

    end do

  end do

  ! penweight=sqrt(penweight)

  call pender(penweight,p,a,lambda,pentype)
  penweight=penweight*real(n,8)

  do count=1,K

    xtemp=x(:,(p*(count-1)+1):(p*count))
    betatemp=beta((p*(count-1)+1):(p*count))
    intvaltemp=intval((m*(count-1)+1):(m*count))
    residualstemp=residuals(:,(m*(count-1)+1):(m*count))

    call cqpen(xtemp,betatemp,intvaltemp,penweight,residualstemp,n,p,tau,m,eps,maxin)

    beta((p*(count-1)+1):(p*count))=betatemp
    intval((m*(count-1)+1):(m*count))=intvaltemp
    residuals(:,(m*(count-1)+1):(m*count))=residualstemp

  end do

  distance=sqrt(sum((beta-beta0)**2)+sum((intval0-intval)**2))

  i=i+1

end do

if((i==maxout) .AND. (distance>eps)) then

  print*,"bad"

end if

return
end subroutine multicqpen



subroutine cqpen(x,beta,intval,penweight,residuals,n,p,tau,m,eps,maxin)
implicit none
integer*4 :: n, p, m, maxin, iter, col, row, nonzerox, i, count3, count4
real*8 :: x(n,p), beta(p), intval(m), penweight(p), residuals(n,m), tau(m), eps, betavecDiff, &
  weight_sum, weight_sum1, weight_sum2, pre_value_final(n*m+1), weight(n*m+1), temp1, temp2, newBeta, &
  betaDiff, pre_int_final(n), weight_int(n)

do iter=1,maxin
  col=0
  betavecDiff=0.d0

  do col=1,p

    weight_sum=0.d0
    nonzerox=0

    do row=1,n

      do i=1,m

        if (x(row,col)/=0.d0) then
          nonzerox=nonzerox+1
          pre_value_final(nonzerox)=residuals(row,i)

          if (pre_value_final(nonzerox)>0.d0) then
            weight(nonzerox)=abs(x(row,col))*(tau(i))
          else
            weight(nonzerox)=abs(x(row,col))*(1-tau(i))
          end if

          pre_value_final(nonzerox)=(pre_value_final(nonzerox)+x(row,col)*beta(col))/(x(row,col))

          weight_sum=weight_sum+weight(nonzerox)
        end if

      end do

    end do

    nonzerox=nonzerox+1
    pre_value_final(nonzerox)=0.d0

    weight(nonzerox)=penweight(col)
    weight_sum=weight_sum+weight(nonzerox)

    call xssort(pre_value_final,weight,nonzerox,2)

    count3=0
    temp1=0.d0
    weight_sum=weight_sum/2.d0
    do while(temp1<=weight_sum )
  
      count3=count3+1
      temp1=temp1+weight(count3)

    end do	

    count4=nonzerox+1
    temp2=0.d0
    do while(temp2<=weight_sum )
  
      count4=count4-1
    
      temp2=temp2+weight(count4)

    end do

    newBeta=(pre_value_final(count3)+pre_value_final(count4))/2.d0

    betaDiff=beta(col)-newBeta
    betavecDiff=betavecDiff+betaDiff*betaDiff
    beta(col)=newBeta

    if (betaDiff/=0.d0) then

      do row=1,n

        do i=1,m

          residuals(row,i)=residuals(row,i)+x(row,col)*betaDiff

        end do

      end do

    end if

  end do

  ! intercepts

    do i=1,m

      weight_sum=0.d0

    do row=1,n

      pre_int_final(row)=residuals(row,i)+intval(i)
      weight_int(row)=1.d0
      weight_sum=weight_sum+weight_int(row)

    end do

    call xssort(pre_int_final,weight_int,n,2)

    count3=0
    temp1=0.d0
    weight_sum1=weight_sum*tau(i)
    do while(temp1<=weight_sum1 )
  
      count3=count3+1
      temp1=temp1+weight_int(count3)

    end do	

    count4=n+1
    temp2=0.d0
    weight_sum2=weight_sum*(1.d0-tau(i))
    do while(temp2<=weight_sum2 )
  
      count4=count4-1
      temp2=temp2+weight_int(count4)

    end do

    newBeta=(pre_int_final(count3)+pre_int_final(count4))/2.d0

    betaDiff=intval(i)-newBeta
    betavecDiff=betavecDiff+betaDiff*betaDiff
    intval(i)=newBeta

    if (betaDiff/=0.d0) then

      do row=1,n

        residuals(row,i)=residuals(row,i)+betaDiff

      end do
 
    end if

  end do

  if(betavecDiff<eps*eps) then

    exit

  end if

end do


return
end subroutine cqpen



Subroutine xssort(x, y, n, kflag)
Integer i, il(21), iu(21), n, nn, kk, kflag, j, m, ij, l, k
Double Precision x(n), y(n), r, t, ty, tty, tt
!***FIRST EXECUTABLE STATEMENT  XSSORT
nn = n
If (nn>=1) Goto 10
Return
10 kk = iabs(kflag)
If ((kk==1) .Or. (kk==2)) Goto 15
Return
!
  ! ALTER ARRAY X TO GET DECREASING ORDER IF NEEDED
!
  15 If (kflag>=1) Goto 30
Do i = 1, nn
x(i) = -x(i)
End Do
30 Goto (100, 200), kk
!
  ! SORT X ONLY
!
  100 Continue
m = 1
i = 1
j = nn
r = .375
110 If (i==j) Goto 155
If (r>.5898437) Goto 120
r = r + 3.90625E-2
Goto 125
120 r = r - .21875
125 k = i
!                                  SELECT A CENTRAL ELEMENT OF THE
!                                  ARRAY AND SAVE IT IN LOCATION T
ij = i + (dfloat(j-i)*r)
t = x(ij)
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER
!                                  THAN T, INTERCHANGE WITH T
If (x(i)<=t) Goto 130
x(ij) = x(i)
x(i) = t
t = x(ij)
130 l = j
!                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
!                                  T, INTERCHANGE WITH T
If (x(j)>=t) Goto 140
x(ij) = x(j)
x(j) = t
t = x(ij)
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER
!                                  THAN T, INTERCHANGE WITH T
If (x(i)<=t) Goto 140
x(ij) = x(i)
x(i) = t
t = x(ij)
Goto 140
135 tt = x(l)
x(l) = x(k)
x(k) = tt
!                                  FIND AN ELEMENT IN THE SECOND HALF OF
!                                  THE ARRAY WHICH IS SMALLER THAN T
140 l = l - 1
If (x(l)>t) Goto 140
!                                  FIND AN ELEMENT IN THE FIRST HALF OF
!                                  THE ARRAY WHICH IS GREATER THAN T
145 k = k + 1
If (x(k)<t) Goto 145
!                                  INTERCHANGE THESE ELEMENTS
If (k<=l) Goto 135
!                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
!                                  THE ARRAY YET TO BE SORTED
If (l-i<=j-k) Goto 150
il(m) = i
iu(m) = l
i = k
m = m + 1
Goto 160
150 il(m) = k
iu(m) = j
j = l
m = m + 1
Goto 160
!                                  BEGIN AGAIN ON ANOTHER PORTION OF
!                                  THE UNSORTED ARRAY
155 m = m - 1
If (m==0) Goto 300
i = il(m)
j = iu(m)
160 If (j-i>=1) Goto 125
If (i==1) Goto 110
i = i - 1
165 i = i + 1
If (i==j) Goto 155
t = x(i+1)
If (x(i)<=t) Goto 165
k = i
170 x(k+1) = x(k)
k = k - 1
If (t<x(k)) Goto 170
x(k+1) = t
Goto 165
!
  ! SORT X AND CARRY Y ALONG
!
  200 Continue
m = 1
i = 1
j = nn
r = .375
210 If (i==j) Goto 255
If (r>.5898437) Goto 220
r = r + 3.90625E-2
Goto 225
220 r = r - .21875
225 k = i
!                                  SELECT A CENTRAL ELEMENT OF THE
!                                  ARRAY AND SAVE IT IN LOCATION T
ij = i + (dfloat(j-i)*r)
t = x(ij)
ty = y(ij)
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER
!                                  THAN T, INTERCHANGE WITH T
If (x(i)<=t) Goto 230
x(ij) = x(i)
x(i) = t
t = x(ij)
y(ij) = y(i)
y(i) = ty
ty = y(ij)
230 l = j
!                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
!                                  T, INTERCHANGE WITH T
If (x(j)>=t) Goto 240
x(ij) = x(j)
x(j) = t
t = x(ij)
y(ij) = y(j)
y(j) = ty
ty = y(ij)
!                                  IF FIRST ELEMENT OF ARRAY IS GREATER
!                                  THAN T, INTERCHANGE WITH T
If (x(i)<=t) Goto 240
x(ij) = x(i)
x(i) = t
t = x(ij)
y(ij) = y(i)
y(i) = ty
ty = y(ij)
Goto 240
235 tt = x(l)
x(l) = x(k)
x(k) = tt
tty = y(l)
y(l) = y(k)
y(k) = tty
!                                  FIND AN ELEMENT IN THE SECOND HALF OF
!                                  THE ARRAY WHICH IS SMALLER THAN T
240 l = l - 1
If (x(l)>t) Goto 240
!                                  FIND AN ELEMENT IN THE FIRST HALF OF
!                                  THE ARRAY WHICH IS GREATER THAN T
245 k = k + 1
If (x(k)<t) Goto 245
!                                  INTERCHANGE THESE ELEMENTS
If (k<=l) Goto 235
!                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
!                                  THE ARRAY YET TO BE SORTED
If (l-i<=j-k) Goto 250
il(m) = i
iu(m) = l
i = k
m = m + 1
Goto 260
250 il(m) = k
iu(m) = j
j = l
m = m + 1
Goto 260
!                                  BEGIN AGAIN ON ANOTHER PORTION OF
!                                  THE UNSORTED ARRAY
255 m = m - 1
If (m==0) Goto 300
i = il(m)
j = iu(m)
260 If (j-i>=1) Goto 225
If (i==1) Goto 210
i = i - 1
265 i = i + 1
If (i==j) Goto 255
t = x(i+1)
ty = y(i+1)
If (x(i)<=t) Goto 265
k = i
270 x(k+1) = x(k)
y(k+1) = y(k)
k = k - 1
If (t<x(k)) Goto 270
x(k+1) = t
y(k+1) = ty
Goto 265
!
  ! CLEAN UP
!
  300 If (kflag>=1) Return
Do i = 1, nn
x(i) = -x(i)
End Do
Return
End Subroutine xssort



subroutine pender(beta,p,a,lambda,pentype)
implicit none
integer*4 :: p, pentype, count
real*8 :: beta(p), a, lambda, temp

! SCAD
do count=1,p
  
  if (pentype==0) then 
		if(abs(beta(count))<lambda) then
				temp=lambda
		else if (abs(beta(count))<a*lambda) then
				temp=(a*lambda-abs(beta(count)))/(a-1.d0)
		else 
				temp=0.d0
		end if
		
	! MCP	
	else if(pentype==1)  then
		
	  if (abs(beta(count))<(a*lambda)) then
		  temp=abs(lambda-(abs(beta(count)/a)))
		else 
			temp=0.d0
		end if
	
	! LASSO	
	else 
	  temp=lambda
	end if
	
	beta(count)=temp
	
end do	

return
end subroutine pender
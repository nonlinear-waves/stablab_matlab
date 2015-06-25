function out=evans(yl,yr,lambda,s,p,m,e)
% out=evans(yl,yr,lambda,s,p,m,e)
%
% Returns the evans function output at a given point.
%
% Input "yl" and "yr" are respectively the initializing values on the left
% and right for the desired manifolds, "lambda" is the value in the complex
% plane where the Evans function is evaluated, and s,p,m,e are structures
% explained in the STABLAB documentation.

fun=str2func(e.evans);
out=fun(yl,yr,lambda,s,p,m,e);

%--------------------------------------------------------------------------
% adj_reg_compound
%--------------------------------------------------------------------------
function out=adj_reg_compound(yl,yr,lambda,s,p,m,e)

Lmani = manifold_compound(e.Li,wedgie(yl),lambda,s,p,m,e.LA,e.kl,1);

Rmani = manifold_compound(e.Ri,wedgie(yr),lambda,s,p,m,e.RA,e.kr,-1);

out=dot(Lmani,Rmani);

%--------------------------------------------------------------------------
% reg_adj_compound
%--------------------------------------------------------------------------
function out=reg_adj_compound(yl,yr,lambda,s,p,m,e)

Lmani = manifold_compound(e.Li,wedgie(yl),lambda,s,p,m,e.LA,e.kl,1);

Rmani = manifold_compound(e.Ri,wedgie(yr),lambda,s,p,m,e.RA,e.kr,-1);

out=dot(Rmani,Lmani);

%--------------------------------------------------------------------------
% reg_adj_polar
%--------------------------------------------------------------------------
function out=reg_adj_polar(WL,WR,lambda,s,p,m,e)

%
% Solve for the basis on left
%

OmegaL0 = orth(WL);
alphaL = OmegaL0'*WL;
muL = trace(OmegaL0'*e.LA(e.Li(1),lambda,s,p)*OmegaL0);
[omegal,gammal] = manifold_polar(e.Li,OmegaL0,lambda,e.LA,s,p,m,e.kl,muL);

%
% Solve for the basis on the right
%

OmegaR0 = orth(WR);
alphaR = OmegaR0'*WR;
muR= trace(OmegaR0'*e.RA(e.Ri(1),lambda,s,p)*OmegaR0);
[omegar,gammar] = manifold_polar(e.Ri,OmegaR0,lambda,e.RA,s,p,m,e.kr,muR);

%
% Evaluate the determinant
%

out = (det(alphaL)*gammal)*conj(det(alphaR)*gammar)*det(omegar'*omegal);

%--------------------------------------------------------------------------
% adj_reg_polar
%--------------------------------------------------------------------------

function out=adj_reg_polar(WL,WR,lambda,s,p,m,e)

%
% Solve for the basis on left
%

OmegaL0 = orth(WL);
alphaL = OmegaL0'*WL;
muL = trace(OmegaL0'*e.LA(e.Li(1),lambda,s,p)*OmegaL0);
[omegal,gammal] = manifold_polar(e.Li,OmegaL0,lambda,e.LA,s,p,m,e.kl,muL);

%
% Solve for the basis on the right
%

OmegaR0 = orth(WR);
alphaR = OmegaR0'*WR;
muR= trace(OmegaR0'*e.RA(e.Ri(1),lambda,s,p)*OmegaR0);
[omegar,gammar] = manifold_polar(e.Ri,OmegaR0,lambda,e.RA,s,p,m,e.kr,muR);

%
% Evaluate the determinant
%

out = conj(det(alphaL)*gammal)*(det(alphaR)*gammar)*det(omegal'*omegar);

%--------------------------------------------------------------------------
% reg_reg_polar
%--------------------------------------------------------------------------

function out = reg_reg_polar(WL,WR,lambda,s,p,m,e)

%
% Solve for the basis on left
%

OmegaL0 = orth(WL);
alphaL = OmegaL0'*WL;
muL = trace(OmegaL0'*e.LA(e.Li(1),lambda,s,p)*OmegaL0);
[omegal,gammal] = manifold_polar(e.Li,OmegaL0,lambda,e.LA,s,p,m,e.kl,muL);

%
% Solve for the basis on the right
%

OmegaR0 = orth(WR);
alphaR = OmegaR0'*WR;
muR= trace(OmegaR0'*e.RA(e.Ri(1),lambda,s,p)*OmegaR0);
[omegar,gammar] = manifold_polar(e.Ri,OmegaR0,lambda,e.RA,s,p,m,e.kr,muR);

%
% Evaluate the determinant
%

out = (det(alphaL)*gammal)*(det(alphaR)*gammar)*det([omegal omegar]);









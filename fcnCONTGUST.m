%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This file is part of jAERO Software
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
%   NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Copyright (C) 2021 by Jan Schwochow (janschwochow@web.de)
%   $Revision: 1.0 $  $Date: 2021/07/19 12:00:00 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cgust,t,Cgust,f] = fcnCONTGUST(V,dt,T,sigma,L)
% calculates random continues vertical gust with
% v.Karman turbulence spectral contents
%
% V     flight speed in m/s
% dt    sample time increment
% T     time length of gust
% cgust continues gust
% t     time vectorize

N = fix(T/dt);
t = [0:1:N-1]*dt;
fs = 1/dt; % sampling frequency
f = [0:1:N-1]/N *fs;
[Om,Phiw] = fcnKARMAN(V,f,sigma,L); % vKarman turbulence spectrum
fw = Om/2/pi*V;
f = fw;
df = (f(2)-f(1))/fs^2;
cgust = fcnFFTNOISE(Phiw'*fs,1);
Cgust = fft(cgust,N); 
cgustrms = sqrt(mean(cgust' .* conj(cgust'))) % root mean square of gust
energy_time = sum(cgust.*conj(cgust) * dt) % energy time domain
energy_freq = sum(Cgust.*conj(Cgust) * df) % energy frequency domain

%     if nargout == 0
        figure(100)
%         h1 = subplot(2,1,1);
        plot(t,real(cgust),'-','linewidth',2);
        grid on, grid minor
        xlabel('time')
        ylabel('gust')
%         h2 = subplot(2,1,2);
%         semilogy(f,abs(Cgust),'linewidth',2);
%         hold on
%         semilogy(fw,abs(Phiw)*fs,'--','linewidth',2);
%         grid on
%         xlim([f(1) 1])
%         ylim([1e-1 max(abs(Cgust))])
%         xlabel('frequency')
%         ylabel('abs gust')
%         legend('gust PSD','vKarman PSD')
%     end
end


function [Om,phiw] = fcnKARMAN(V,f,sigma,L)
% power spectral density for the v.Karman turbulence spectrum
% L    = 25;   % turbulence scale in m
%L    = 2500;  % turbulence scale in ft
Om   = 2*pi*f/V; % reduced frequency in rad/m
phiw = (sigma*sigma)*L/pi*(1+8/3*(1.339*Om*L).^2)./...
            (1+    (1.339*Om*L).^2).^(11/6);  % von Karman spectrum
end


function noise = fcnFFTNOISE(f,Nseries)
% Generate noise with a given power spectrum.
% f         the fft of a time series (must be a column vector)
% Nseries   number of noise series to generate. (default=1)
% noise     surrogate series with same power spectrum as f.
if nargin<2; Nseries = 1; end
f = f(:);
N = length(f);
Np = floor((N-1)/2);
phases = rand(Np,Nseries)*2*pi;
phases = complex(cos(phases),sin(phases));
f = repmat(f,1,Nseries);
f(2:Np+1,:) = f(2:Np+1,:).*phases;
f(end:-1:end-Np+1,:) = conj(f(2:Np+1,:));
noise = real(ifft(f,[],1));
end

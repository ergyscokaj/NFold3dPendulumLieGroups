clc;
clear all;
close all;

prompt = 'How many 3D pendulums do you want to connect?\n\n';
P = 2; %input(prompt);

L = rand(P,1); 
m = rand(P,1);

L = 0*L + 1; 
m = 0*m + 1;

[q0,w0,z0] = initializeSE3N(P);
% [q0,w0,z0] = initializeSE3N_smallVariation(P);
% [q0,w0,z0] = initializeSE3N_largeVariation(P);

Energy = @(q,w) 0.5*w'*assembleR(q,L,m)*w + potential(q,L,m);

disp("Energy of this initial condition: "+num2str(Energy(q0,w0)));

t0 = 0;
T = 5; 
N = 1000; 
time = linspace(t0,T,N); 
dt = time(2)-time(1);

getq = @(v) extractq(v);
getw = @(v) extractw(v);

f = @(v) fManiToAlgebra(getq(v),getw(v),L,m); 
action = @(B,input) actionSE3N(B,input); 
vecField = @(sigma,p) dexpinvSE3N(sigma,f(action(exponentialSE3N(sigma),p)));


z = z0;
qC = zeros(3*P,N);
pC = qC;

Len = zeros(3*P,1);
for i = 1:P
    Len(3*i-2:3*i) = L(i)*ones(3,1);
end
Mat = diag(Len);
if P>1
    for i = 3:3:3*(P-1)
        Mat = Mat + diag(Len(1:3*P-i),-i);
    end
end
qC(:,1) = q0;
pC(:,1) = Mat*q0;

%% TIME EVOLUTION OF THE SOLUTION

prompt = 'Do you want to see the Time Evolution of the solution? Write 1 for yes, 0 for no\n\n';
C1 = 0;%input(prompt);

zC = zeros(6*P,N);
if C1==1
    zC(:,1) = z0;
    for i = 1:N-1
        zC(:,i+1) = FreeRK4SE3N(f,action,dt,zC(:,i));

        qC(:,i+1) = extractq(zC(:,i+1));
        pC(:,i+1) = Mat*qC(:,i+1);
    end

    figure('Units','normalized','Position',[0 0 1 1])
    t = 0;
    for i = 1:2:N
        t = dt*i;
        plot3([0,pC(1,i)],[0,pC(2,i)],[0,pC(3,i)],'k-*',...
            [pC(3*(1:P-1)-2,i),pC(3*(1:P-1)+1,i)],[pC(3*(1:P-1)-1,i),pC(3*(1:P-1)+2,i)],...
                [pC(3*(1:P-1),i),pC(3*(1:P-1)+3,i)],'k-*','Markersize',5, 'LineWidth',3);
        xlabel("x")
        ylabel("y")
        zlabel("z")
        str = "Time evolution of the pendulums, T="+num2str(t);
        axis([-sum(L) sum(L) -sum(L) sum(L) -sum(L) sum(L)]);
        title(str)
        pause(0.00000000000001);
    end
end

    %% PRESERVATION OF THE GEOMETRY

if P==2
    prompt = 'Do you want to see the Time Evolution of norms of the solution? Write 1 for yes, 0 for no\n\n';
    C2 = 1;%input(prompt);

    if C2==1

        [nn,NN] = compareNorms(f,action,vecField,z0,L,m);
        times = linspace(0,5,NN);
        
%         % Set the default font to 'Times New Roman'
%          set(0,'DefaultTextFontName', 'Times New Roman')
%          set(0,'DefaultAxesFontName', 'Times New Roman')
%         % Set the default font size for math expressions
%          set(0,'DefaultTextInterpreter','latex')

         subplot(4, 2, 1)
         plot(times, 1-nn(:,:,1).^2,'-','linewidth',2);
%          yticks(1)
         xlim([0 5])
         legend('$1-q_1^\top q_1$','$1-q_2^\top q_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Norm",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('ODE45','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 2)
         plot(times, 1-nn(:,:,2).^2,'-','linewidth',2);
         xlim([0 5])
%          yticks(1)
         legend('$1-q_1^\top q_1$','$1-q_2^\top q_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Norm",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('RK4','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 3)
         plot(times, 1-nn(:,:,3).^2,'-','linewidth',2);
         xlim([0 5])
%          yticks(1)
         legend('$1-q_1^\top q_1$','$1-q_2^\top q_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Norm",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('Lie Euler','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 4)
         plot(times, 1-nn(:,:,4).^2,'-','linewidth',2);
         xlim([0 5])
%          yticks(1)
         legend('$1-q_1^\top q_1$','$1-q_2^\top q_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Norm",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('Lie Euler Heun','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 5)
         plot(times, 1-nn(:,:,5).^2,'-','linewidth',2);
         xlim([0 5])
%          yticks(1)
         legend('$1-q_1^\top q_1$','$1-q_2^\top q_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Norm",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('CF4','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 6)
         plot(times, 1-nn(:,:,6).^2,'-','linewidth',2);
         xlim([0 5])
%          yticks(1)
         legend('$1-q_1^\top q_1$','$1-q_2^\top q_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Norm",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('RKMK4','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 7)
         plot(times, 1-nn(:,:,7).^2,'-','linewidth',2);
         xlim([0 5])
%          yticks(1)
         legend('$1-q_1^\top q_1$','$1-q_2^\top q_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Norm",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('RKMK4 with 2 commutators','fontsize',25);
         set(ax,'TickLength',[0 0]);

         subplot(4, 2, 8)
         plot(times, 1-nn(:,:,8).^2,'-','linewidth',2);
         xlim([0 5])
%          yticks(1)
         legend('$1-q_1^\top q_1$','$1-q_2^\top q_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Norm",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('RKMK3','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
        
         
%          h = sgtitle("$1-q_i(t)^\top q_i(t)$",'interpreter','latex');
        sgtitle("$1-q_i(t)^\top q_i(t)$","FontSize",30,'interpreter','latex');
%         set(gcf, 'Position',  [0  50 500 700]);
    end




    prompt = 'Do you want to see the Time Evolution of the scalar products w'' * q? Write 1 for yes, 0 for no\n\n';
    C3 = 1;%input(prompt);

    if C3==1

        [nn,NN] = tangentBehaviour(f,action,vecField,z0,L,m);
        times = linspace(0,5,NN);
        
        n1 = nn(:,:,1);
        n2 = nn(:,:,2);
        n3 = nn(:,:,3);
        n4 = nn(:,:,4);
        n5 = nn(:,:,5);
        n6 = nn(:,:,6);
        n7 = nn(:,:,7);
        n8 = nn(:,:,8);
        
        figure;
         subplot(4, 2, 1)
         plot(times, nn(:,:,1),'-','linewidth',2);
         xlim([0 5])
         legend('$q_1^\top\omega_1$','$q_2^\top \omega_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Product",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('ODE45','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 2)
         plot(times, nn(:,:,2),'-','linewidth',2);
         xlim([0 5])
         legend('$q_1^\top\omega_1$','$q_2^\top \omega_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Product",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('RK4','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 3)
         plot(times, nn(:,:,3),'-','linewidth',2);
         xlim([0 5])
         legend('$q_1^\top\omega_1$','$q_2^\top \omega_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Product",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"]) 
         title('Lie Euler','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 4)
         plot(times, nn(:,:,4),'-','linewidth',2);
         xlim([0 5])
         legend('$q_1^\top\omega_1$','$q_2^\top \omega_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Product",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('Lie Euler Heun','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 5)
         plot(times, nn(:,:,5),'-','linewidth',2);
         xlim([0 5])
         legend('$q_1^\top\omega_1$','$q_2^\top \omega_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Product",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('CF4','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 6)
         plot(times, nn(:,:,6),'-','linewidth',2);
         xlim([0 5])
         legend('$q_1^\top\omega_1$','$q_2^\top \omega_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Product",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('RKMK4','fontsize',25);
         set(ax,'TickLength',[0 0]);
         
         subplot(4, 2, 7)
         plot(times, nn(:,:,7),'-','linewidth',2);
         xlim([0 5])
         legend('$q_1^\top\omega_1$','$q_2^\top \omega_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25);
         ylabel("Product",'fontsize',25);
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('RKMK4 with 2 commutators','fontsize',25);
         set(ax,'TickLength',[0 0]);

         subplot(4, 2, 8)
         plot(times, nn(:,:,8),'-','linewidth',2);
         xlim([0 5])
         legend('$q_1^\top\omega_1$','$q_2^\top \omega_2$','interpreter','latex','FontSize',20,'Location', 'southeast', Color = "#eaeaf3");
         xlabel("Time",'fontsize',25,'Color', 'k');
         ylabel("Product",'fontsize',25,'Color', 'k');
         ax = gca;
         ax.Color = "#eaeaf3";
         ax.XAxis.FontSize = 25;
         ax.YAxis.FontSize = 25;
         ax.XGrid = 'on';
         ax.YGrid = 'on';
         ax.GridColorMode = 'manual';
         ax.GridColor = 'w';
         ax.GridAlpha = 1;
         ax.GridLineStyle = '-'; 
         ax.YColor = 'k';
         ax.XColor = 'k';

         colororder(["#1f77b4";"#ff7f0e";"#2ca02c";"#d62728";"#9467bd";"#8c564b";"#e377c2";"#7f7f7f";"#bcbd22";"#17becf"])
         title('RKMK3','fontsize',25);

         sgtitle("$q_i(t)^\top\omega_i(t)$","FontSize",30,'interpreter','latex');
         set(ax,'TickLength',[0 0]);
    end
end

%% CONVERGENCE RATE OF THE METHODS

prompt = 'Do you want to see the convergence rate of the implemented methods? Write 1 for yes, 0 for no\n\n';
C = 1;%input(prompt);
if C==1
    checkConvergenceRate(f,action,vecField,z0,L,m); 
end



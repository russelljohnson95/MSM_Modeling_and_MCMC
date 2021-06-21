function pairs_v3(x,panelfun,names,skip,varargin)
%PAIRS pairs plot like in R
% PAIRS(X,PANELFUN,NAMES) pairs plot of matrix X
% PANELFUN(X(:,i),X(:,j)) if exists is applied to every panel
% NAMES is a shell array of column names

% ML 2000

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:38 $
[n,p] = size(x);

lb = [0.5, 0, 0, 4, 0, 0.05, -0.2];
ub = [2, 0.5, 6, 18, 0.06, 0.15, 0.2];

if p>50
   error('too many pairs')
end

if nargin<4 | isempty(skip)
  skip=1;
end
inds=1:skip:n;
count = 1;
%clf
for j=2:p
  for i=1:j-1
    if p==2
      h=gca;
    else
      h=subplot(p-1,p-1,(j-2)*(p-1)+i);
    end
    plot(x(inds,i),x(inds,j),'.');
    [R,P(:,:,j)] = corrcoef(x(inds,i),x(inds,j));
    txt = ['r = ',num2str(R(2,1),2)]; 
    x1min = min(x(inds,i))-abs(min(x(inds,i))*0.2);
    x1max = max(x(inds,i))+abs(max(x(inds,i))*0.1); 
    x2min = min(x(inds,j))-abs(min(x(inds,j))*0.2); 
    x2max = max(x(inds,j))+abs(max(x(inds,j))*0.1); 
    xlim([lb(i),ub(i)])
    ylim([lb(j),ub(j)])
    yticks([lb(j) ub(j)])
    yticklabels([lb(j) ub(j)])
    xticks([lb(i) ub(i)])
    xticklabels([lb(i) ub(i)])

%     xlim([(x1min),(x1max)])
%     ylim([(x2min),(x2max)])
%     yticks([floor(x2min) ceil(x2max)])
%     yticklabels([floor(x2min) ceil(x2max)])
%     xticks([floor(x1min) ceil(x1max)])
%     xticklabels([floor(x1min) ceil(x1max)])
    x1(i) = x1min;
    x2(j) = x2max*1.15; 
%     text(x1(i),x2(j),txt,'fontsize',13,'Color','red');
    title(txt)
    set(gca,'fontsize',8)
%     pbaspect([3 2 1])
    if j~=p
      set(h,'xtick',[])      
    end
    if i~=1
      set(h,'ytick',[])
      
    end
    box off
%     ytickangle(45)
    if i==1 & nargin>2 & ~isempty(names)
      ylabel(names{j},'FontSize',10,'FontWeight','bold')
    end
    if j == p 
        xlabel(names{i},'FontSize',10','FontWeight','bold')
    end
%      set(gca, 'YTick', []);
     
     pvalue(count) = P(2,1);
     count = count+1; 
%      ylim([-20 20])
%      xlim([-20 20])
%     if i==j-1 & nargin>2 & ~isempty(names)
%       if p==2
%         xlabel(names{i},'FontSize',16);
%       else
%         title(names{i},'FontSize',16)
%       end
%     end
%     if nargin>1 & ~isempty(panelfun)
%       feval(panelfun,x(:,i),x(:,j),varargin{:});
%     end
  end   
end

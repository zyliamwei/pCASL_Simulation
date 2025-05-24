function y=ZLap(data,range,step,delta)
% return the laplace domain information
% data: original to be Laplace transformed
% range: range of Laplace domain
% step: step of the Laplace domain
% delta: dwelling time of the data
if nargin~=4
    error('Invalid input ...');
else
    if min(size(data))~=1
        error('Invalid matrix dimension ...');
    else
        RD=min(range);
        RU=max(range);
        Seq=RD:step:RU;
        Scale=length(Seq);
        SeqData=0:length(data)-1;
        ReData=zeros(1,Scale);
        for ni=1:length(Seq)
            ReData(ni)=sum(data.*exp(Seq(ni)*SeqData*delta));
        end
    end
end
y=ReData;
end
            
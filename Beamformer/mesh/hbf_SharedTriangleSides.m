function [stri1,stri0,stri_more]=hbf_SharedTriangleSides(input)
%HBF_SHAREDTRIANGLESIDES finds the connections between sides of triangles.
%
% function [stri1,stri0,stri_more]=SHAREDTRIANGLESIDES(mesh)
% function [stri1,stri0,stri_more]=SHAREDTRIANGLESIDES(elements)
%   mesh:   triangle mesh, hbf struct
%   elements:   triangle description, [N(triangles) x 3]
%
%   stri1.tri: triangles that share sides, [N(triangles) x 3]
%   stri1.sside: which side of the target triangle is shared, [N(triangles) x 3]%
%   stri0: list of triangles that have a side that is not shared
%   stri_more: list of triangles that have a side that is shared with more
%     than one other triangle
%
%   If there are more than two triangles that contain the same side, stri1
%   has -1 in the corresponding fields
%   If there is a side that belongs to only one triangle, stri1 and sside1
%   have value of 0 for that side.
%   This BEM solver needs a closed, well-behaving mesh. That is, each side 
%   of each triangle is shared with exactly one side of another triangle.
%   In other words, stri0 and stri_more should be empty.
%
%   Sides are organized as follows: for triangle p1, p2, p3,
%       side1=p2-p1, side2=p3-p1, side3=p3-p2
%
% v160229 Matti Stenroos
if isstruct(input)
    e=input.e;
    noe=size(e,1);
else
    e=input;
    noe=size(e,1);
end
stri1=zeros(noe,3);
sside1=zeros(noe,3);
s1=sort(e(:,[1 2]),2);
s2=sort(e(:,[1 3]),2);
s3=sort(e(:,[2 3]),2);

%make unique IDs for alle sides, so comparison is quicker
mcoef=noe;
s1=mcoef*s1(:,1)+s1(:,2);
s2=mcoef*s2(:,1)+s2(:,2);
s3=mcoef*s3(:,1)+s3(:,2);
s=[s1 s2 s3];
tridone=zeros(noe,3);


% stri_more.sourcetri=[];stri_more.sourceside=[];
% stri_more.targettri={};stri_more.targetside={};
% stri0.tri=[];
% stri0.side=[];
stri_more=[];
stri0=[];
% Nsofar=0;
for E=1:noe,
    for SR=1:3,
        if tridone(E,SR)
            continue;
        end
        ref=s(E,SR);
        tester=any(s==ref,2);
        tester(E)=0;
        if any(tester)
            targettri=find(tester);
            Ntri=size(targettri,1);
            if Ntri==1,
                targetside=find(s(tester,:)==ref);
                stri1(E,SR)=targettri;
                sside1(E,SR)=targetside;
                stri1(targettri,targetside)=E;
                sside1(targettri,targetside)=SR;
                tridone(E,SR)=1;
                tridone(targettri,targetside)=1;
           
            else
                stri_more=[stri_more;E];
                stri1(E,SR)=-1;
                sside1(E,SR)=-1;
%                 targetside=zeros(Ntri,1);
%                 for I=1:Ntri,
%                     targetside(I)=find(s(targettri(I),:)==ref);
%                 end
%                 stri_more.sourcetri=[stri_more.sourcetri;E];
%                 stri_more.sourceside=[stri_more.sourceside;SR];
%                 Nsofar=length(stri_more.sourceside);
%                 stri_more.targettri{Nsofar}=[E targettri'];
%                 stri_more.targetside{Nsofar}=[SR targetside'];
                tridone(E,SR)=1;
            end
        else
            stri1(E,SR)=0;
            sside1(E,SR)=0;
%             stri0.tri=[stri0.tri;E];
%             stri0.side=[stri0.side;SR];
            stri0=[stri0;E];
            tridone(E,SR)=1;
            
        end
    end
end
% if Nsofar==0,
%     stri_more=[];
% end
% if isempty(stri0.tri)
%     stri0=[];
% end

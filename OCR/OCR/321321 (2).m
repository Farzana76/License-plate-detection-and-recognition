
image=imread('image.jpg')
for a=0:1000
    for b=0:1000
        for c=0:1000
            if (image(a,b,c)<=40)
                image(a,b,c)=0;
            else
            end
        end
    end
end
image
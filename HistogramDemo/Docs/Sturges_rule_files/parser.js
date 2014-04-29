// JavaScript Document

function parser(fn) {
var X, Y, a,b,g,tt, link;

a = location.href.lastIndexOf('\\');
b = document.location.href;
var pat = "ModelRiskHelp";
Y = b.lastIndexOf(pat);
Y=Y+pat.length;
tt = b.substring(0,Y);
X = fn.lastIndexOf('Models');
if (X==-1)
{ 
  X = fn.lastIndexOf('Videos');
};
g = fn.substring(0,6)+'/'+ fn.substring(6,fn.length);
g = fn;
/////
if (b.lastIndexOf(".chm")!=-1 && fn.lastIndexOf('Videos')!=-1)
{
   g = g.substring(0,g.length-4)+".exe";
   //alert(g);
}
/////
fn = g;
link = document.location.href.substring(0, Y+1) + fn;
if (b.lastIndexOf(".chm")!=-1)
{
  var X1, Y1, sl, a1, ra;
  ra = /:/;
  a1 = location.href.search(ra);
  if (a1 == 2)
  X1 = 14;
  else
  X1 = 7;
  sl = "\\";
  Y1 = location.href.lastIndexOf(sl) + 1;
  link = 'file:///' + location.href.substring(X1, Y1) + fn;
  //location.href = link;  
 
  //alert(link); 
  //link = link.substring(link.lastIndexOf(":\\")-1,link.length);
}
if (fn.lastIndexOf(".swf")!=-1)
{
 //alert(link);  
 window.open(link, '_blank', 'width=800,height=600');
}
else
{
   if (fn.lastIndexOf(".htm")!=-1)
   {
      //alert(link);  
      window.open(link, '_blank', 'width=800,height=600');
   }
   else
   { 
    //alert(link);
    location.href = link;
   }
}
}


// script
function frame_loader(y)
{

var frm = window.parent.frames;
    
    var index_page = "index.htm#";
    var pat = "ModelRiskHelp";
    var k = document.location.href.lastIndexOf(".chm");
    if(frm.length==0 && k==-1)
    {
        var pos1 = document.location.href.lastIndexOf(pat);
        var st1 = document.location.href.substring(0,pos1+pat.length+1);
        var st2 = document.location.href.substring(pos1+pat.length+1,document.location.href.length); 
        var targetURL=st1+index_page + st2;
        //alert(targetURL);
   
        document.location.href = targetURL;
    }
}


// JavaScript Document

//function parser(fn) {
//var X, Y, sl, a, ra, link;
//ra = /:/;
//a = location.href.search(ra);
//if (a == 2)
//X = 14;
//else
//X = 7;
//sl = "\\";
//Y = location.pathname.lastIndexOf(sl) + 1;
//link = 'file:///' + location.href.substring(X, Y) + fn;
//location.href = link;
//}

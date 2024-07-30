$(window).bind("pageshow", function() {
    $(".userAnswer").removeAttr("disabled");
    document.getElementById('userAnswer').value = "";
    $("#togglertext").click();
});

var exons = [];

$(document).ready(function(){
    $('img').click(function() {
        $("[id^=radio]").prop("checked",false);
        var id = "#".concat(this.id)
        var val = parseInt(this.id)
        var temp = ''
        exons.push(val)
        //$(id).fadeTo(1000,0.1);
        var largest = 0;
        var smallest = 80;
        for (i=0; i<=exons.length;i++){
            if (exons[i]>largest) {
                var largest = parseInt(exons[i]);
            }
            if (exons[i]<smallest) {
                smallest = parseInt(exons[i]);
            }
        }
        if (smallest === largest){
            temp = "Deletion of exon "
            temp = temp.concat(this.id)
        }
        else if (smallest !== largest){
            temp = "Deletion of exons "
            temp = temp +smallest + "-"+ largest
        }
        for (i=smallest;i<=largest;i++){
            id = "#" + i
            $(id).fadeTo(1000,0.1);
        }
        $(".userAnswer").val(temp)
        $(".userAnswer").attr("disabled","disabled");
    });

    $("#reset").click(function(){
        $("[id^=radio]").prop("checked",false);
        $("img").fadeTo(1000,1);
        $(".userAnswer").val("");
        $(".userAnswer").removeAttr("disabled");
        exons.length = 0;
    });

    $(".btn-primary").click(function(){
        $(".userAnswer").removeAttr("disabled");
    });
    $("[id^=radio]").click(function(){
        $(".userAnswer").val($(this).val());
    });
    $(".userAnswer").click(function(){
        $("[id^=radio]").prop("checked",false);
    });
});

function handleUserAnswer(e){
    var el = e.target || e.srcElement;
    document.getElementById('.userAnswer').value = document.getElementById('.userAnswer').value.concat(el.value);
}

var radios = document.getElementsByTagName('input');
for(var i = 0; i<radios.length; i++){
    var r = radios[i];
    if(r.getAttribute('name') == 'userAnswerChoice'){
        if(r.addEventListener){
            r.addEventListener('change',handleUserAnswer);
        }else if(r.attachEvent){
            r.attachEvent('onchange',handleUserAnswer);
        }else{
            r.onChange = handleUserAnswer;
        }
    }
}


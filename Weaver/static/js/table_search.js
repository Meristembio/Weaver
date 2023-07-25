function do_filter_default(){
    do_filter($("#table_search-input").first().val());
}
function do_filter(value){
    value = value.replace(/[^A-Za-z0-9]/g,'').toLowerCase();
    $(".table-search-target tbody tr").filter(function() {
        var element_value = $(this).find('.table-search-search_on').first().attr('data-search').replace(/[^A-Za-z0-9]/g,'').toLowerCase();
        console.log(element_value + " / " + value);
        if(element_value.indexOf(value) > -1){
            $(this).removeClass('table_search-hide');
        } else {
            $(this).addClass('table_search-hide');
        }
    });
}

$(document).ready(function(){
    $("#table_search-input").on("keyup paste", function(e) {
        var element = $(this);
        var event = e;
        setTimeout(function() {
            if (event.which == 27) { // Esc
                element.val("");
            }
            do_filter(element.val());
        }, 100);
    });
    $("#table_search-clear").click(function() {
        $(".table-search-target tbody tr").removeClass('table_search-hide');
        $("#table_search-input").val('');
    });
});
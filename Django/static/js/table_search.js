var search_on = 'all';

function do_filter_default(){
    do_filter($("#table_search-input").first().val());
}

function do_filter(value){
    value = value.replace(/[^A-Za-z0-9]/g,'').toLowerCase();
    if(value){
        $(".table-search-target tbody tr").filter(function() {
            var element = $(this).find('.table-search-search_on').first().attr('data-search-' + search_on)
            if(search_on != 'idx'){
                var file_names = $(this).find('[data-sanger-file-name]').map(function(){
                    return $(this).attr('data-sanger-file-name') || '';
                }).get().join(' ');
                element = (element || '') + ' ' + file_names;
            }
            if(element){
                var element_value = element.replace(/[^A-Za-z0-9]/g,'').toLowerCase();
                if(search_on == 'idx'){
                    // Perfect match
                    if(element_value == value){
                        $(this).removeClass('table_search-hide');
                    } else {
                        $(this).addClass('table_search-hide');
                    }
                } else {
                    // Contains
                    if(element_value.indexOf(value) > -1){
                        $(this).removeClass('table_search-hide');
                    } else {
                        $(this).addClass('table_search-hide');
                    }
                }
            }
        });
    }
    else {
        // SHow all
        $(".table-search-target tbody tr").removeClass('table_search-hide');
    }
}

$(document).ready(function(){
    if ($('[data-server-search="true"]').length) {
        return;
    }
    $('.table_search-target').click(function(){
        $('.table_search-target').removeClass('active');
        $(this).addClass('active');
        search_on = $(this).attr('data-type');
        do_filter($('#table_search-input').val());
    });
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
